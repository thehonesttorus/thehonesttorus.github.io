"""Goal inference as correspondence synthesis (Section 6).

Goals are second-quantized predicates of words from the library (C6, fock-ladder): for a
grade-k predicate p, dGamma(p) counts the k-tuples satisfying p and Gamma(p) asserts all do.
Families (Section 6.1), ranked by the conditional Parry kernel at beta_c:

    grade 0:  Zero[c]  count_c = 0                        d = 16
              Eq[c,c'] count_c = count_c' (both seen)       d = 120
    grade 1:  Absent[t] / Exist[t]  dGamma(1_t) = 0 / >= 1  d = |Xi_0| each
    grade 2:  Some[R,t1,t2]  dGamma(p_R) >= 1 ; All[R,t1,t2] Gamma(p_R)   d = 6 |Xi_0|^2 each
              R in {contact, coincide, adjacent1}

Evidence (6.2): every observed non-winning post-state is a certified negative (a goal true
there dies); the post-state of a completing action is never observed, so positives are
evaluated on the model's predicted post-state when the action's coherent set agrees, else
kept pending. Before any positive, the set is the prior restricted to words false on every
frame observed so far, ranked by omega(g) * A(g) with A(g) the fraction of g's observables
that any stored transition has moved (6.3, Component 4). Goal-splitting (6.4) targets the
reachable state maximizing W(s)(W_tot - W(s)).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Optional

import numpy as np

from .algebra import BaseAlgebra, Component
from .library import Library
from .store import Store

RELATIONS = ("contact", "coincide", "adjacent1")


def _relation(R: str, a: Component, b: Component, A: BaseAlgebra) -> bool:
    if R == "contact":
        return (min(a.id, b.id), max(a.id, b.id)) in A.contact
    if R == "coincide":
        return a.bbox == b.bbox and a.g == b.g
    if R == "adjacent1":
        return abs(a.bbox[0] - b.bbox[0]) + abs(a.bbox[1] - b.bbox[1]) == 1 and a.g == b.g
    raise ValueError(R)


@dataclass(frozen=True)
class GoalWord:
    family: str
    params: tuple

    def eval(self, A: BaseAlgebra) -> bool:
        f, p = self.family, self.params
        if f == "Zero":
            return int(A.hist[p[0]]) == 0
        if f == "Eq":
            c1 = sum(1 for c in A.comps if c.color == p[0]); c2 = sum(1 for c in A.comps if c.color == p[1])
            return c1 == c2
        if f == "Exist":
            return any(c.key == p[0] for c in A.comps)
        if f == "Absent":
            return not any(c.key == p[0] for c in A.comps)
        if f in ("Some", "All"):
            R, t1, t2 = p
            xs = [c for c in A.comps if c.key == t1]
            ys = [c for c in A.comps if c.key == t2]
            if not xs or not ys:
                return False
            hits = [any(_relation(R, x, y, A) for y in ys if y.id != x.id) for x in xs]
            return any(hits) if f == "Some" else all(hits)
        raise ValueError(f)

    def dim(self, xi0: int) -> float:
        return {"Zero": 16.0, "Eq": 120.0, "Exist": float(xi0), "Absent": float(xi0),
                "Some": 6.0 * xi0 * xi0, "All": 6.0 * xi0 * xi0}[self.family]

    def observables(self) -> set:
        f, p = self.family, self.params
        if f == "Zero":
            return {("colour", p[0])}
        if f == "Eq":
            return {("colour", p[0]), ("colour", p[1])}
        if f in ("Exist", "Absent"):
            return {("type", p[0])}
        return {("type", p[1]), ("type", p[2]), ("rel", p[0], p[1], p[2])}

    def describe(self) -> str:
        return f"{self.family}{self.params if self.family in ('Zero','Eq') else '[types]'}"


@dataclass
class GoalSet:
    lib: Library
    G: dict[GoalWord, float] = field(default_factory=dict)      # goal -> log weight (prior)
    negatives_level: list[BaseAlgebra] = field(default_factory=list)
    positives: list[BaseAlgebra] = field(default_factory=list)
    pending: list[tuple] = field(default_factory=list)            # (A_T, a_T, data_T)
    confirmed: int = 0
    killed: int = 0
    regenerated: int = 0

    # ---- vocabulary --------------------------------------------------------------
    def _alphabet(self, frames: Iterable[BaseAlgebra]) -> tuple[list[int], list[bytes]]:
        colours, keys = set(), set()
        counts: dict[bytes, int] = {}
        for A in frames:
            colours.update(int(c) for c in np.nonzero(A.hist)[0].tolist())
            for c in A.frame_vectors():
                counts[c.key] = counts.get(c.key, 0) + 1
        keys = sorted(counts, key=lambda k: -counts[k])
        return sorted(colours), keys

    def _prior(self, frames: list[BaseAlgebra]) -> dict[GoalWord, float]:
        colours, keys = self._alphabet(frames)
        xi0 = max(2, len(keys))
        words: list[GoalWord] = []
        words += [GoalWord("Zero", (c,)) for c in colours]
        words += [GoalWord("Eq", (c1, c2)) for i, c1 in enumerate(colours) for c2 in colours[i + 1:]]
        words += [GoalWord("Exist", (k,)) for k in keys] + [GoalWord("Absent", (k,)) for k in keys]
        pair_keys = keys[:12]
        if pair_keys:
            for R in RELATIONS:
                for t1 in pair_keys:
                    for t2 in pair_keys:
                        words.append(GoalWord("Some", (R, t1, t2)))
                        words.append(GoalWord("All", (R, t1, t2)))
        return {w: -self.lib.beta_c * float(np.log(w.dim(xi0))) for w in words}

    # ---- evidence ----------------------------------------------------------------
    def negative(self, A: BaseAlgebra) -> None:
        self.negatives_level.append(A)
        if len(self.negatives_level) > 400:
            self.negatives_level = self.negatives_level[-400:]
        if self.G:
            dead = [g for g in self.G if g.eval(A)]
            for g in dead:
                del self.G[g]
            self.killed += len(dead)

    def positive(self, A_T: BaseAlgebra, a_T: int, data_T: Optional[dict]) -> None:
        self.pending.append((A_T, a_T, data_T))

    def resolve_pending(self, H: dict, predicted_frame_fn) -> None:
        still = []
        for (A_T, a_T, data_T) in self.pending:
            hyps = H.get(a_T) or []
            if not hyps:
                still.append((A_T, a_T, data_T)); continue
            frames = set()
            posts = []
            for h in hyps[:4]:
                fr = predicted_frame_fn(h.cs, A_T, data_T)
                if fr is None:
                    continue
                k = np.ascontiguousarray(fr, dtype=np.int8).tobytes()
                if k not in frames:
                    frames.add(k); posts.append(fr)
            if len(posts) != 1:
                still.append((A_T, a_T, data_T)); continue
            from .algebra import build
            A_post = build(posts[0])
            self.positives.append(A_post)
            self.confirmed += 1
            if self.G:
                self.G = {g: w for g, w in self.G.items() if g.eval(A_post)}
        self.pending = still

    def new_level(self) -> None:
        self.negatives_level = []

    def ensure(self, A_now: BaseAlgebra) -> None:
        """Regenerate from the prior when empty: consistent with all positives and this level's negatives."""
        if self.G:
            return
        frames = self.negatives_level[-50:] + self.positives + [A_now]
        prior = self._prior(frames)
        G = {}
        for g, w in prior.items():
            if any(g.eval(P) is False for P in self.positives):
                continue
            if g.eval(A_now):
                continue
            if any(g.eval(N) for N in self.negatives_level[-200:]):
                continue
            G[g] = w
        self.G = G
        self.regenerated += 1

    # ---- ranking and splitting -----------------------------------------------------
    def alignment(self, store: Store) -> dict[GoalWord, float]:
        moved: set = set()
        for recs in store.records.values():
            for r in recs:
                for i, j, g in r.tr.moved:
                    moved.add(("type", r.A0.comps[i].key)); moved.add(("colour", r.A0.comps[i].color))
                for i, j in r.tr.recolours:
                    moved.add(("colour", r.A0.comps[i].color)); moved.add(("colour", r.A1.comps[j].color))
                    moved.add(("type", r.A0.comps[i].key)); moved.add(("type", r.A1.comps[j].key))
                for i in r.tr.annihilated:
                    moved.add(("type", r.A0.comps[i].key)); moved.add(("colour", r.A0.comps[i].color))
                for j in r.tr.created:
                    moved.add(("type", r.A1.comps[j].key)); moved.add(("colour", r.A1.comps[j].color))
                if r.A0.contact != r.A1.contact:
                    moved.add("contact")
        out = {}
        for g in self.G:
            obs = g.observables()
            hit = 0
            for o in obs:
                if o[0] == "rel":
                    if "contact" in moved or ("type", o[2]) in moved or ("type", o[3]) in moved:
                        hit += 1
                elif o in moved:
                    hit += 1
            out[g] = hit / max(len(obs), 1)
        return out

    def weights(self, store: Optional[Store] = None) -> dict[GoalWord, float]:
        if not self.G:
            return {}
        logs = np.array(list(self.G.values()))
        w = np.exp(logs - logs.max())
        W = dict(zip(self.G.keys(), w))
        if self.positives or store is None:
            return W
        al = self.alignment(store)
        return {g: W[g] * (0.5 + 0.5 * al[g]) for g in W}   # alignment ranks, it does not kill (6.3)

    def split_target(self, reach: dict, store: Optional[Store] = None) -> Optional[list]:
        """Path to the reachable state maximizing (W(W_tot - W), W, -len)."""
        W = self.weights(store)
        if not W:
            return None
        W_tot = sum(W.values())
        best = None
        for key, (path, A) in reach.items():
            if not path:
                continue
            Ws = sum(w for g, w in W.items() if g.eval(A))
            if Ws <= 0:
                continue
            k = (Ws * (W_tot - Ws), Ws, -len(path))
            if best is None or k > best[0]:
                best = (k, path)
        return best[1] if best else None
