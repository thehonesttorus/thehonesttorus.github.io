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


def _index(A: BaseAlgebra) -> tuple[np.ndarray, dict]:
    """Per-algebra index for goal evaluation: component counts per colour, components per shape key."""
    idx = A.__dict__.get("_gidx")
    if idx is None:
        cc = np.bincount([c.color for c in A.comps], minlength=16) if A.comps else np.zeros(16, dtype=np.int64)
        bk: dict = {}
        for c in A.comps:
            bk.setdefault(c.key, []).append(c)
        idx = (cc, bk)
        A.__dict__["_gidx"] = idx
    return idx


def _facts(A: BaseAlgebra) -> dict:
    """Every goal word's truth value on A, as sets (cached per algebra): colours with mass 0, colour pairs with equal
    component counts, shape types present, and for each relation R the type pairs (t1, t2) with some related pair and
    with every t1 component related to some t2 component. Goal words evaluate by lookup, and a level's negatives
    are summarised by unions of these sets, so no negative is ever dropped."""
    f = A.__dict__.get("_gfacts")
    if f is not None:
        return f
    cc, bk = _index(A)
    zero = frozenset(c for c in range(16) if int(A.hist[c]) == 0)
    eq = frozenset((c1, c2) for c1 in range(16) for c2 in range(c1 + 1, 16) if int(cc[c1]) == int(cc[c2]))
    keys = frozenset(bk)
    partners: dict[tuple, set] = {}
    def rel(R: str, a: Component, b: Component) -> None:
        partners.setdefault((R, a.key, b.key), set()).add(a.id)
        partners.setdefault((R, b.key, a.key), set()).add(b.id)
    for i, j in A.contact:
        rel("contact", A.comps[i], A.comps[j])
    by_pose: dict[tuple, list] = {}
    for c in A.comps:
        by_pose.setdefault((c.bbox, c.g), []).append(c)
    for group in by_pose.values():
        for x in range(len(group)):
            for y in range(x + 1, len(group)):
                rel("coincide", group[x], group[y])
    by_corner: dict[tuple, list] = {}
    for c in A.comps:
        by_corner.setdefault((c.bbox[0], c.bbox[1], c.g), []).append(c)
    for c in A.comps:
        for dy, dx in ((1, 0), (0, 1)):
            for o in by_corner.get((c.bbox[0] + dy, c.bbox[1] + dx, c.g), []):
                rel("adjacent1", c, o)
    some = frozenset(partners)
    allr = frozenset(k for k, ids in partners.items() if len(ids) == len(bk.get(k[1], ())))
    f = dict(zero=zero, eq=eq, keys=keys, some=some, all=allr)
    A.__dict__["_gfacts"] = f
    return f


@dataclass(frozen=True)
class GoalWord:
    family: str
    params: tuple

    def eval(self, A: BaseAlgebra) -> bool:
        f, p = self.family, self.params
        if f == "Zero":
            return int(A.hist[p[0]]) == 0
        cc, bk = _index(A)
        if f == "Eq":
            return int(cc[p[0]]) == int(cc[p[1]])
        if f == "Exist":
            return p[0] in bk
        if f == "Absent":
            return p[0] not in bk
        if f in ("Some", "All"):
            facts = _facts(A)
            return tuple(p) in (facts["some"] if f == "Some" else facts["all"])
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
    negatives_level: list[BaseAlgebra] = field(default_factory=list)   # the latest few, for the alphabet only
    neg_union: dict = field(default_factory=dict)                     # unions of goal facts over every negative
    key_counts: dict = field(default_factory=dict)                    # shape-type frequencies over the level's negatives
    colours_seen: set = field(default_factory=set)
    positives: list[BaseAlgebra] = field(default_factory=list)
    pending: list[tuple] = field(default_factory=list)            # (A_T, a_T, data_T)
    confirmed: int = 0
    killed: int = 0
    regenerated: int = 0

    # ---- vocabulary --------------------------------------------------------------
    def _alphabet(self, frames: Iterable[BaseAlgebra]) -> tuple[list[int], list[bytes]]:
        colours, keys = set(self.colours_seen), set()
        counts: dict[bytes, int] = dict(self.key_counts)
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
        pair_keys = keys
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
        if len(self.negatives_level) > 8:
            self.negatives_level = self.negatives_level[-8:]
        f = _facts(A)
        u = self.neg_union
        if not u:
            u.update(zero=set(), eq=set(), keys_any=set(), keys_all=set(f["keys"]), some=set(), all=set())
        u["zero"] |= f["zero"]; u["eq"] |= f["eq"]; u["keys_any"] |= f["keys"]; u["keys_all"] &= f["keys"]
        u["some"] |= f["some"]; u["all"] |= f["all"]
        for c in A.frame_vectors():
            self.key_counts[c.key] = self.key_counts.get(c.key, 0) + 1
        self.colours_seen.update(int(c) for c in np.nonzero(A.hist)[0].tolist())
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
            for h in hyps:                         # the whole coherent set must agree (spec 6.2)
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
        self.neg_union = {}
        self.key_counts = {}
        self.colours_seen = set()

    def _killed_by_negatives(self, g: "GoalWord") -> bool:
        """True iff g is true on some negative of this level (read from the unions; exact for every family)."""
        u = self.neg_union
        if not u:
            return False
        f, p = g.family, g.params
        if f == "Zero":
            return p[0] in u["zero"]
        if f == "Eq":
            return tuple(sorted(p)) in u["eq"]
        if f == "Exist":
            return p[0] in u["keys_any"]
        if f == "Absent":
            return p[0] not in u["keys_all"]
        if f == "Some":
            return tuple(p) in u["some"]
        if f == "All":
            return tuple(p) in u["all"]
        return False

    def ensure(self, A_now: BaseAlgebra) -> None:
        """Regenerate from the prior when empty: consistent with all positives and this level's negatives."""
        if self.G:
            return
        frames = self.negatives_level + self.positives + [A_now]
        prior = self._prior(frames)
        G = {}
        for g, w in prior.items():
            if any(g.eval(P) is False for P in self.positives):
                continue
            if g.eval(A_now):
                continue
            if self._killed_by_negatives(g):
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
