"""Planning inside the verified model (Section 7; N16: never over the true environment).

Uniform-cost best-first search over predicted frames, using one consensus hypothesis per
action id (the top-weight member of its coherent set; disagreement is resolved by the Kato
rule, not here). Death is priced by replay length, never forbidden. The visited set is keyed
by the predicted frame (the only hash use besides the ignorance table, N2). Successor base
algebras are built lazily through a cache, so the per-step cost is bounded by the time
budget rather than by the branching factor. `plan_reachable` returns every state reached
with its path, for goal-splitting (Section 6.4).
"""
from __future__ import annotations

import heapq
import time
from typing import Callable, Optional

import numpy as np

from . import algebra, kato
from .words import Hypothesis, predicted

Action = tuple[int, Optional[dict]]


def consensus(H: dict[int, list[Hypothesis]]) -> dict[int, Hypothesis]:
    return {a: hs[0] for a, hs in H.items() if hs}


def _key(fr: np.ndarray) -> bytes:
    return np.ascontiguousarray(fr, dtype=np.int8).tobytes()


def plan_reachable(A0: algebra.BaseAlgebra, avail: list[int], Hc: dict[int, Hypothesis], node_budget: int = 60,
                   time_budget_s: float = 0.08, perceive=None, memory=None) -> dict[bytes, tuple[list[Action], algebra.BaseAlgebra]]:
    """Reachable predicted states within the budget: frame key -> (path, base algebra)."""
    t0 = time.time()
    start_key = A0.key()
    reach: dict[bytes, tuple[list[Action], algebra.BaseAlgebra]] = {start_key: ([], A0)}
    frontier = [(0.0, 0, start_key)]
    counter = 1
    expanded = 0
    while frontier and expanded < node_budget and time.time() - t0 < time_budget_s:
        cost, _, key = heapq.heappop(frontier)
        path, A = reach[key]
        expanded += 1
        for (a, data) in kato.candidates(A, avail):
            h = Hc.get(a)
            if h is None or h.cs.over:
                continue
            fr, flag = predicted(h.cs, A, data, memory)
            if fr is None or flag == "GAME_OVER":
                continue
            k2 = _key(fr)
            if k2 in reach:
                continue
            A2 = perceive(fr) if perceive is not None else algebra.build_cached(fr)
            reach[k2] = (path + [(a, data)], A2)
            counter += 1
            heapq.heappush(frontier, (cost + 1.0, counter, k2))
            if time.time() - t0 > time_budget_s:
                break
    return reach


def plan_to(A0: algebra.BaseAlgebra, avail: list[int], Hc: dict[int, Hypothesis], target: Callable[[algebra.BaseAlgebra], bool],
            node_budget: int = 60, time_budget_s: float = 0.08) -> Optional[list[Action]]:
    reach = plan_reachable(A0, avail, Hc, node_budget, time_budget_s)
    best: Optional[list[Action]] = None
    for key, (path, A) in reach.items():
        if path and target(A) and (best is None or len(path) < len(best)):
            best = path
    return best
