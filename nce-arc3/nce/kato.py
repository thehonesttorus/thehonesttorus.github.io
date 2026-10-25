"""The Kato exploration rule (Section 5): KMS-expected tangent tax of the coherent set.

D(a, data) = sum_{h<h'} p_h p_h' sin^2 theta(y_h, y_h') = 1/2 (1 - Tr rho^2), with p_h the
normalized KMS weights of the hypotheses for (a, data) and cos theta = agreeing cells / 4096
between their predicted frames (Theorem 4.1(i) read on rank-one prediction projections; the
E10 commutator form 1/2 sin^2 2theta is its small-angle equivalent and has a saddle, so the
tax is used). Identical predictions merge into one class.

Order of business: untried actions first (every outcome newborn, T2); then argmax of
D / cost over candidates, cost = 1 + P_death * replay (C8); an action with no hypotheses
scores only in a newborn context (its (type multiset, pointed type) never stored), with
D = 1/2 (1 - 1/m) over the m applicable letter families (moment-map span rule); a
non-deterministic action scores zero. Returns None when every candidate has D = 0.
"""
from __future__ import annotations

from typing import Optional

import numpy as np

from .algebra import BaseAlgebra
from .store import Store
from .words import ClauseSet, Hypothesis, predict, predicted, predicted_frame, word_context

KATO_BUDGET_S = 0.25      # compute bound per decision (seconds); candidates are visited most diverse first
def families_applicable(lib=None) -> int:
    """m in the newborn deficit 1/2 (1 - 1/m): the letter families that can act on a selected component, i.e. the
    library's edges leaving the component vertex C, plus the identity. Read from the library graph, not set by hand."""
    if lib is None:
        from .library import Library
        lib = Library(win_levels=2)
    return 1 + sum(1 for p in lib.graph.primitives.values() if p.source == "C")


FAMILIES_APPLICABLE = families_applicable()


def candidates(A: BaseAlgebra, avail: list[int]) -> list[tuple[int, Optional[dict]]]:
    out: list[tuple[int, Optional[dict]]] = [(a, None) for a in avail if a not in (0, 6)]
    if 6 in avail:
        # one representative cell per component, every component: first one per shape type (HN order), then the rest
        comps = sorted(A.comps, key=lambda c: (c.stratum, c.id))
        seen: set = set()
        first, rest = [], []
        for c in comps:
            (rest if c.type in seen else first).append(c)
            seen.add(c.type)
        for c in first + rest:
            x, y = c.representative_cell()
            out.append((6, {"x": x, "y": y}))
    return out


def context_key(A: BaseAlgebra, data: Optional[dict]) -> tuple:
    """The context of an experiment before any word exists: the multiset of frame-vector types and, for a pointed
    action, the pointed component's type and pose (its gauge element: bbox corner and D4 element). Two instances of
    one type differ by a gauge element, and whether the outcome depends on it is exactly what is untested."""
    types = tuple(sorted((k, n) for k, n in A.type_counts_fv().items()))
    pointed = None
    if data is not None:
        cid = int(A.cell2comp[int(data["y"]), int(data["x"])])
        if cid >= 0:
            c = A.comps[cid]
            pointed = (c.color, c.key, c.offset, c.g)
    return (types, pointed)


def record_context(r) -> tuple:
    """context_key of a stored demonstration (cached on the record; cleared on re-perception)."""
    k = r.__dict__.get("_ctx")
    if k is None:
        k = context_key(r.A0, r.data)
        r.__dict__["_ctx"] = k
    return k


def newborn_context(a: int, A: BaseAlgebra, data: Optional[dict], store: Store, seen_cache: dict) -> bool:
    key = context_key(A, data)
    seen = seen_cache.get(a)
    if seen is None:
        seen = {record_context(r) for r in store.transitions(a)}
        seen_cache[a] = seen
    return key not in seen


def representative_data(a: int, A: BaseAlgebra) -> Optional[dict]:
    if a != 6 or not A.comps:
        return None
    # first component in HN order of the lowest-slope strata (frame vectors) when any exists, else the top stratum
    non_top = [c for c in A.comps if c.stratum > 0]
    c = sorted(non_top or A.comps, key=lambda c: (-c.stratum, c.id))[0]
    x, y = c.representative_cell()
    return {"x": x, "y": y}


def purity_deficit(hyps: list[Hypothesis], A: BaseAlgebra, data: Optional[dict], memory=None) -> float:
    """1/2 (1 - Tr rho^2) over the classes of distinct predicted frames, weighted by KMS weights; the overlap of two
    classes is the flat one-hot metric, cos = (cells where the predicted frames agree) / 4096."""
    classes: dict[bytes, list] = {}
    logws = np.array([h.logw for h in hyps])
    w = np.exp(logws - logws.max())
    for h, wi in zip(hyps, w):
        fr, flag = predicted(h.cs, A, data, memory)
        if fr is None:
            key = b"GAME_OVER" if flag == "GAME_OVER" else b"OFFGRID"
        else:
            fr = fr.ravel()
            key = np.ascontiguousarray(fr, dtype=np.int8).tobytes()
        if key in classes:
            classes[key][0] += wi
        else:
            classes[key] = [wi, fr]
    if len(classes) < 2:
        return 0.0
    items = list(classes.values())
    Z = sum(c[0] for c in items)
    D = 0.0
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            fi, fj = items[i][1], items[j][1]
            if fi is None or fj is None:
                cos2 = 0.0
            else:
                cos = 1.0 - float(np.count_nonzero(fi != fj)) / (64 * 64)
                cos2 = cos * cos
            D += (items[i][0] / Z) * (items[j][0] / Z) * (1.0 - cos2)
    return float(D)


def kato_select(A: BaseAlgebra, avail: list[int], H: dict[int, list[Hypothesis]], partials: dict[int, list[ClauseSet]],
                store: Store, actions_since_level: int, restrict: Optional[set] = None, memory=None) -> Optional[tuple[int, Optional[dict]]]:
    # 1. try once: an action with no demonstration (RESET is learned like any other, but only when offered)
    for a in avail:
        if a == 0:
            continue
        if not store.tried(a):
            return (a, representative_data(a, A))
    seen_cache: dict = {}
    best, bestval, besttie = None, 0.0, (0.0, 0)
    import time as _time
    t_start = _time.time()
    for (a, data) in candidates(A, avail):
        if best is not None and _time.time() - t_start > KATO_BUDGET_S:
            break
        if restrict is not None and (a, None if data is None else (data["x"], data["y"])) not in restrict:
            continue
        if a in store.nondeterministic:
            continue
        post = store.seen_post(a, A, data)
        if post == b"GAME_OVER":
            continue          # a certified death in exactly this context teaches nothing (C8)
        if post is not None and post == A.key():
            continue          # a certified no-op in exactly this context teaches nothing (C8)
        hyps = H.get(a) or []
        from_partials = False
        if not hyps and partials.get(a):
            hyps = [Hypothesis(cs, 0.0) for cs in partials[a]]
            from_partials = True
        if hyps and not from_partials:
            # newborn relative to the words: some member has never been verified on what it sees here
            newborn = any(word_context(h.cs, A, data) not in h.contexts for h in hyps)
        else:
            newborn = newborn_context(a, A, data, store, seen_cache)
        if hyps and not from_partials:
            D = purity_deficit(hyps, A, data, memory)
            p_death = sum(np.exp(h.logw) for h in hyps if h.cs.over or predict(h.cs, A, data)[2] == "GAME_OVER") / max(sum(np.exp(h.logw) for h in hyps), 1e-300)
            if newborn:
                # the word has never been tested in this context: its unfixed branch is a direct sum (T7)
                D = max(D, 0.5 * (1.0 - 1.0 / FAMILIES_APPLICABLE))
        else:
            if not newborn:
                continue
            D = 0.5 * (1.0 - 1.0 / FAMILIES_APPLICABLE)
            p_death = 0.0
        cost = 1.0 + p_death * max(actions_since_level, 1)
        val = D / cost
        # ties (equal commutator value) go to the least-demonstrated action: the ignorance table's emptiest row
        if val <= 0.0:
            continue
        tie = (val, -len(store.transitions(a)))
        if best is None or tie > besttie:
            best, bestval, besttie = (a, data), val, tie
    return best
