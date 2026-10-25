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
from .words import ClauseSet, Hypothesis, predicted_frame

CANDIDATE_CAP = 64
FAMILIES_APPLICABLE = 6   # identity, Step/Move, Rot, Recolour, Annihilate, Create (Fill unfixed inside each)


def candidates(A: BaseAlgebra, avail: list[int]) -> list[tuple[int, Optional[dict]]]:
    out: list[tuple[int, Optional[dict]]] = [(a, None) for a in avail if a not in (0, 6)]
    if 6 in avail:
        # one representative cell per component; HN order (top stratum first), capped
        comps = sorted(A.comps, key=lambda c: (c.stratum, c.id))[:CANDIDATE_CAP]
        for c in comps:
            x, y = c.representative_cell()
            out.append((6, {"x": x, "y": y}))
    return out


def context_key(A: BaseAlgebra, data: Optional[dict]) -> tuple:
    types = tuple(sorted((k, n) for k, n in A.type_counts().items()))
    pointed = None
    if data is not None:
        cid = int(A.cell2comp[int(data["y"]), int(data["x"])])
        pointed = A.comps[cid].key if cid >= 0 else None
    return (types, pointed)


def newborn_context(a: int, A: BaseAlgebra, data: Optional[dict], store: Store, seen_cache: dict) -> bool:
    key = context_key(A, data)
    seen = seen_cache.get(a)
    if seen is None:
        seen = {context_key(r.A0, r.data) for r in store.transitions(a)}
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


def purity_deficit(hyps: list[Hypothesis], A: BaseAlgebra, data: Optional[dict]) -> float:
    """1/2 (1 - Tr rho^2) over the classes of distinct predicted frames, weighted by KMS weights."""
    classes: dict[bytes, list] = {}
    logws = np.array([h.logw for h in hyps])
    w = np.exp(logws - logws.max())
    for h, wi in zip(hyps, w):
        fr = predicted_frame(h.cs, A, data)
        if fr is None:
            key = b"OFFGRID"
            diff = None
        else:
            diff = fr.ravel() != A.frame.ravel()
            key = np.packbits(diff).tobytes()
        if key in classes:
            classes[key][0] += wi
        else:
            classes[key] = [wi, diff]
    if len(classes) < 2:
        return 0.0
    items = list(classes.values())
    Z = sum(c[0] for c in items)
    D = 0.0
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            di, dj = items[i][1], items[j][1]
            if di is None or dj is None:
                cos2 = 0.0
            else:
                cos = 1.0 - float(np.count_nonzero(di != dj)) / (64 * 64)
                cos2 = cos * cos
            D += (items[i][0] / Z) * (items[j][0] / Z) * (1.0 - cos2)
    return float(D)


def kato_select(A: BaseAlgebra, avail: list[int], H: dict[int, list[Hypothesis]], partials: dict[int, list[ClauseSet]],
                store: Store, actions_since_level: int, restrict: Optional[set] = None) -> Optional[tuple[int, Optional[dict]]]:
    # 1. try once: an action with no demonstration (RESET is learned like any other, but only when offered)
    for a in avail:
        if a == 0:
            continue
        if not store.tried(a):
            return (a, representative_data(a, A))
    seen_cache: dict = {}
    best, bestval = None, 0.0
    for (a, data) in candidates(A, avail):
        if restrict is not None and (a, None if data is None else (data["x"], data["y"])) not in restrict:
            continue
        if a in store.nondeterministic:
            continue
        hyps = H.get(a) or []
        if not hyps and partials.get(a):
            hyps = [Hypothesis(cs, 0.0) for cs in partials[a]]
        if hyps:
            D = purity_deficit(hyps, A, data)
            p_death = sum(np.exp(h.logw) for h in hyps if h.cs.over) / max(sum(np.exp(h.logw) for h in hyps), 1e-300)
        else:
            if not newborn_context(a, A, data, store, seen_cache):
                continue
            D = 0.5 * (1.0 - 1.0 / FAMILIES_APPLICABLE)
            p_death = 0.0
        cost = 1.0 + p_death * max(actions_since_level, 1)
        val = D / cost
        if val > bestval:
            best, bestval = (a, data), val
    return best
