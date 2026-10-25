"""Decomposition of an observed transition into gauge moves, recolourings, creations and annihilations.

Theory: Proposition 1.3 (gauge torus action: the transition splits into charge sectors by
colour fibre), eq. (1.3) (creation / annihilation on the Fock ladder), Theorem 1.5(iii)
(the KMS weight of a move is set by the statistical dimension of its family), T8 (the same
thing moved = identity composed with a gauge element), C11 (transform, then match: first
moments give the translation, the canonical form gives the D4 element, fibres give the colour).

Given two base algebras A0 (before) and A1 (after):
  charge[c]   = count of colour-c components after minus before
  moves       = (i, j, Gauge) pairs with equal shape type: a translation by (dy, dx) and a D4
                change; weight set by the Parry prior of the translation family (4*l members on
                the l1-sphere of radius l, l >= 1) and of D4 (7 non-identity elements)
  recolours   = (i, j) pairs with equal shape key, different colour, and overlapping cells
  annihilated = components of A0 explained by nothing
  created     = components of A1 explained by nothing

The assignment inside a type bucket maximises total KMS weight (minimises summed cost), so
"prefer small displacements" is the Parry order, not a heuristic.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from .algebra import BaseAlgebra, Component, GRID

try:
    from scipy.optimize import linear_sum_assignment as _lsa  # type: ignore
except Exception:  # pragma: no cover
    _lsa = None


@dataclass(frozen=True)
class Gauge:
    """A gauge element of D4 x Z^2 relating a component before to a component after."""
    dy: int
    dx: int
    g_before: int
    g_after: int

    @property
    def l1(self) -> int:
        return abs(self.dy) + abs(self.dx)

    @property
    def is_identity(self) -> bool:
        return self.dy == 0 and self.dx == 0 and self.g_before == self.g_after

    @property
    def rotated(self) -> bool:
        return self.g_before != self.g_after


STEP_DIM, MOVE_DIM, ROT_DIM = 4.0, 127.0 * 127.0 - 1.0, 7.0   # library dimensions of Step, Move, Rot


def translation_cost(l1, beta: float = 1.0):
    """-log Parry weight of the cheapest word for a translation of l1-length l: Step^l (4^l summands) or Move (16128)."""
    l1 = np.asarray(l1, dtype=np.float64)
    return beta * np.where(l1 > 0, np.minimum(l1 * np.log(STEP_DIM), np.log(MOVE_DIM)), 0.0)


def move_cost(gauge: Gauge, beta: float = 1.0) -> float:
    """-log of the Parry weight of the word explaining the gauge element: identity 0; a translation by its cheapest
    word (Step^l or Move); a D4 change beta log 7; both: the serial composition."""
    c = float(translation_cost(gauge.l1, beta)) if gauge.l1 > 0 else 0.0
    if gauge.rotated:
        c += beta * np.log(ROT_DIM)
    return float(c)


@dataclass
class Transition:
    charge: np.ndarray                                  # (16,) component-count change per colour
    moves: list[tuple[int, int, Gauge]] = field(default_factory=list)      # (id0, id1, gauge)
    recolours: list[tuple[int, int]] = field(default_factory=list)         # (id0, id1)
    reshaped: list[tuple[int, int]] = field(default_factory=list)          # (id0, id1) same colour, overlapping
    inplace_rot: list[tuple[int, int, Gauge]] = field(default_factory=list)  # reshaped pairs that are also congruent under a D4 change
    occluded: list[int] = field(default_factory=list)                      # ids in A0 fully covered by moved/created components in A1
    annihilated: list[int] = field(default_factory=list)                   # ids in A0
    created: list[int] = field(default_factory=list)                       # ids in A1
    changed_cells: int = 0                               # Hamming distance between the frames

    @property
    def moved(self) -> list[tuple[int, int, Gauge]]:
        return [(i, j, g) for i, j, g in self.moves if not g.is_identity]

    @property
    def is_null(self) -> bool:
        return self.changed_cells == 0


def _assign(cost: np.ndarray) -> list[tuple[int, int]]:
    if cost.size == 0:
        return []
    if _lsa is not None:
        r, c = _lsa(cost)
        return list(zip(r.tolist(), c.tolist()))
    # greedy fallback: repeatedly take the cheapest remaining pair
    pairs: list[tuple[int, int]] = []
    used_r: set[int] = set()
    used_c: set[int] = set()
    order = np.argsort(cost, axis=None)
    for flat in order.tolist():
        i, j = divmod(flat, cost.shape[1])
        if i in used_r or j in used_c:
            continue
        pairs.append((i, j))
        used_r.add(i)
        used_c.add(j)
        if len(used_r) == cost.shape[0] or len(used_c) == cost.shape[1]:
            break
    return pairs


def _containment_shift(ma: np.ndarray, a: Component, b: Component) -> Optional[tuple[int, int]]:
    """Translation v with transported(a) subset of b, via FFT cross-correlation of the masks (None if none)."""
    mb = b.mask()
    fa = np.fft.rfft2(ma.astype(np.float32))
    fb = np.fft.rfft2(mb.astype(np.float32))
    corr = np.fft.irfft2(fb * np.conj(fa), s=ma.shape)
    # candidate shifts where the overlap equals |a| (full containment), tolerance for float error
    cand = np.argwhere(corr > a.n - 0.5)
    if cand.size == 0:
        return None
    best = None
    for dy, dx in cand.tolist():
        dy = dy - 64 if dy > 32 else dy
        dx = dx - 64 if dx > 32 else dx
        y0, x0 = a.bbox[0] + dy, a.bbox[1] + dx
        if y0 < 0 or x0 < 0 or a.bbox[2] + dy >= GRID or a.bbox[3] + dx >= GRID:
            continue
        shifted = np.zeros_like(ma)
        shifted[y0:a.bbox[2] + dy + 1, x0:a.bbox[3] + dx + 1] = ma[a.bbox[0]:a.bbox[2] + 1, a.bbox[1]:a.bbox[3] + 1]
        if np.all(mb[shifted]):
            l1 = abs(dy) + abs(dx)
            if best is None or l1 < best[0]:
                best = (l1, (dy, dx))
    return best[1] if best else None


def _gauge_between(a: Component, b: Component) -> Gauge:
    return Gauge(b.bbox[0] - a.bbox[0], b.bbox[1] - a.bbox[1], a.g, b.g)


def decompose(A0: BaseAlgebra, A1: BaseAlgebra, beta: float = 1.0, max_shift: int = GRID,
              prefer: Optional[frozenset] = None) -> Transition:
    """Match components of A1 to components of A0 and classify the remainder."""
    charge = np.bincount([c.color for c in A1.comps], minlength=16) - np.bincount([c.color for c in A0.comps], minlength=16)
    tr = Transition(charge=charge.astype(np.int64), changed_cells=int(np.count_nonzero(A0.frame != A1.frame)))
    matched0: set[int] = set()
    matched1: set[int] = set()

    # 1. same shape type (colour, canonical key). Candidate matchings of each bucket are charged by the word that
    #    explains them: each distinct gauge element once (one clause SelType . body covers every pair it displaces),
    #    the identity and the displacements this action has already shown for the type (a clause already in its word)
    #    free. Three candidates: the per-pair assignment, history-first, shared-first; the cheapest charge wins.
    prefer = prefer or frozenset()
    for t, ids0 in A0.by_type.items():
        ids1 = A1.by_type.get(t)
        if not ids1:
            continue
        c0 = [A0.comps[i] for i in ids0]
        c1 = [A1.comps[j] for j in ids1]
        bb0 = np.array([c.bbox for c in c0], dtype=np.int64)
        bb1 = np.array([c.bbox for c in c1], dtype=np.int64)
        g0 = np.array([c.g for c in c0], dtype=np.int64)
        g1 = np.array([c.g for c in c1], dtype=np.int64)
        dy = bb1[None, :, 0] - bb0[:, None, 0]
        dx = bb1[None, :, 1] - bb0[:, None, 1]
        l1 = np.abs(dy) + np.abs(dx)
        rot = g1[None, :] != g0[:, None]
        cost = translation_cost(l1, beta) + np.where(rot, beta * np.log(ROT_DIM), 0.0) + 1e-6 * l1   # tie: the nearer
        cost[l1 > max_shift] = 1e9
        inplace: set[tuple[int, int]] = set()
        touch = rot & (l1 <= max_shift) & (bb0[:, None, 0] <= bb1[None, :, 2]) & (bb1[None, :, 0] <= bb0[:, None, 2]) \
            & (bb0[:, None, 1] <= bb1[None, :, 3]) & (bb1[None, :, 1] <= bb0[:, None, 3])
        for r, c in zip(*np.nonzero(touch)):
            if np.intersect1d(c0[r].cells, c1[c].cells, assume_unique=True).size > 0:
                # a congruent component overlapping its old cells under a D4 change: the cheapest word is the
                # identity with a changed rendering (cost 0: occlusion, reveal, a trail), so it is paired here as
                # a reshape; the Rot reading (cost log 7) is kept in `inplace_rot` for synthesis to try next.
                cost[r, c] = 0.0
                inplace.add((int(r), int(c)))
        known = np.zeros_like(rot)
        if prefer:
            for r in range(len(c0)):
                for c in range(len(c1)):
                    if not rot[r, c] and (t, int(dy[r, c]), int(dx[r, c])) in prefer:
                        known[r, c] = True
        cost_k = np.where(known, 1e-6 * l1, cost)

        def complete(fixed: list) -> list:
            fr = {r for r, _ in fixed}
            fc = {c for _, c in fixed}
            cm = cost_k.copy()
            if fr:
                cm[list(fr), :] = 1e9
                cm[:, list(fc)] = 1e9
            out = list(fixed)
            for r, c in _assign(cm):
                if cm[r, c] < 1e9 and r not in fr and c not in fc:
                    out.append((int(r), int(c)))
            return out

        def charge(pairs: list) -> float:
            elems: dict = {}
            for r, c in pairs:
                if (r, c) in inplace or known[r, c] or (dy[r, c] == 0 and dx[r, c] == 0 and not rot[r, c]):
                    continue
                key = (int(dy[r, c]), int(dx[r, c]), bool(rot[r, c]))
                elems[key] = float(cost[r, c])
            return sum(elems.values()) + 1e-6 * float(sum(l1[r, c] for r, c in pairs))

        def by_vector(mask) -> list:
            rr, cc = np.nonzero(mask & ~rot & (l1 <= max_shift))
            vs: dict = {}
            for r_, c_ in zip(rr.tolist(), cc.tolist()):
                vs.setdefault((int(dy[r_, c_]), int(dx[r_, c_])), []).append((r_, c_))
            fixed, ur, uc = [], set(), set()
            for v in sorted(vs, key=lambda v: (-len(vs[v]), float(translation_cost(abs(v[0]) + abs(v[1]), beta)))):
                for r_, c_ in vs[v]:
                    if r_ not in ur and c_ not in uc:
                        fixed.append((r_, c_)); ur.add(r_); uc.add(c_)
            return fixed

        candidates = [complete([])]
        if known.any():
            candidates.append(complete(by_vector(known)))
        if len(c0) > 1 and len(c1) > 1:
            shared_mask = np.zeros_like(rot)
            rr, cc = np.nonzero(~rot & (l1 <= max_shift))
            counts: dict = {}
            for r_, c_ in zip(rr.tolist(), cc.tolist()):
                counts[(int(dy[r_, c_]), int(dx[r_, c_]))] = counts.get((int(dy[r_, c_]), int(dx[r_, c_])), 0) + 1
            for r_, c_ in zip(rr.tolist(), cc.tolist()):
                if counts[(int(dy[r_, c_]), int(dx[r_, c_]))] >= 2:
                    shared_mask[r_, c_] = True
            if shared_mask.any():
                candidates.append(complete(by_vector(shared_mask)))
        best_pairs = min(candidates, key=charge)
        for r, c in best_pairs:
            gg = Gauge(int(dy[r, c]), int(dx[r, c]), int(g0[r]), int(g1[c]))
            if (r, c) in inplace:
                tr.reshaped.append((ids0[r], ids1[c]))
                tr.inplace_rot.append((ids0[r], ids1[c], gg))
            else:
                tr.moves.append((ids0[r], ids1[c], gg))
            matched0.add(ids0[r])
            matched1.add(ids1[c])

    # 2. recolourings: same canonical key, different colour, overlapping cells (exact position)
    left0 = [c for c in A0.comps if c.id not in matched0]
    left1 = [c for c in A1.comps if c.id not in matched1]
    for a in left0:
        best: Optional[Component] = None
        for b in left1:
            if b.id in matched1 or b.key != a.key or b.color == a.color:
                continue
            if b.bbox == a.bbox and b.g == a.g:
                best = b
                break
        if best is not None:
            tr.recolours.append((a.id, best.id))
            matched0.add(a.id)
            matched1.add(best.id)

    # 2b. containment moves: an unmatched component whose cells, transported by a translation v,
    #     lie inside a same-colour component of A1 (it merged with a neighbour or was partly occluded).
    #     v is read from the cross-correlation of the two masks (second-order statistics, C11), then
    #     confirmed exactly by containment. D4 changes are not searched here.
    left0 = [c for c in A0.comps if c.id not in matched0 and c.stratum > 0]
    left1 = [c for c in A1.comps if c.id not in matched1]
    if left0 and left1:
        for a in left0:
            best = None
            ma = a.mask()
            for b in left1:
                if b.id in matched1 or b.color != a.color or b.n < a.n:
                    continue
                v = _containment_shift(ma, a, b)
                if v is None:
                    continue
                gg = Gauge(v[0], v[1], a.g, a.g)
                if gg.l1 == 0 or gg.l1 > max_shift:
                    continue
                c = move_cost(gg, beta)
                if best is None or c < best[0]:
                    best = (c, b, gg)
            if best is not None:
                tr.moves.append((a.id, best[1].id, best[2]))
                matched0.add(a.id)
                # the target component may absorb several movers; mark it matched only once all is painted
                matched1.add(best[1].id) if best[1].n == a.n else None

    # 3. reshaped: same colour, overlapping cells, different canonical key (e.g. the top stratum,
    #    whose shape is the complement of everything else, or a component that grew or split).
    #    Matched by maximal overlap; the correspondence is annihilation+creation at the library
    #    level unless a word explains it (the matcher only records the pairing).
    left0 = [c for c in A0.comps if c.id not in matched0]
    left1 = [c for c in A1.comps if c.id not in matched1]
    if left0 and left1:
        m1 = np.full(GRID * GRID, -1, dtype=np.int64)
        for b in left1:
            m1[b.cells] = b.id
        cand: list[tuple[int, int, int]] = []
        for a in left0:
            hits = m1[a.cells]
            hits = hits[hits >= 0]
            if hits.size == 0:
                continue
            ids, counts = np.unique(hits, return_counts=True)
            for bid, cnt in zip(ids.tolist(), counts.tolist()):
                if A1.comps[bid].color == a.color:
                    cand.append((int(cnt), a.id, bid))
        cand.sort(reverse=True)
        for cnt, aid, bid in cand:
            if aid in matched0 or bid in matched1:
                continue
            tr.reshaped.append((aid, bid))
            matched0.add(aid)
            matched1.add(bid)

    # 4. occlusion: an unmatched A0 component whose cells are all covered in A1 by components that moved there
    #    or appeared is hidden, not annihilated (strata persist under occlusion; the frame is their rendering).
    unmatched0 = [c for c in A0.comps if c.id not in matched0]
    unmatched1 = [c.id for c in A1.comps if c.id not in matched1]
    if unmatched0:
        cover = np.zeros(GRID * GRID, dtype=bool)
        for i, j, g in tr.moves:
            if not g.is_identity:
                cover[A1.comps[j].cells] = True
        top = len(A1.strata) > 1
        for i, j in tr.reshaped:
            if top and A1.comps[j].stratum == 0:
                continue          # the top stratum's complement always flows into vacated cells; it hides nothing
            cover[A1.comps[j].cells] = True
        for j in unmatched1:
            cover[A1.comps[j].cells] = True
        for c in unmatched0:
            if np.all(cover[c.cells]):
                tr.occluded.append(c.id)
                matched0.add(c.id)
    tr.annihilated = [c.id for c in A0.comps if c.id not in matched0]
    tr.created = unmatched1
    return tr


def summarize(tr: Transition, A0: BaseAlgebra, A1: BaseAlgebra) -> str:
    parts = []
    for i, j, g in tr.moved:
        parts.append(f"move c{A0.comps[i].color} n{A0.comps[i].n} by ({g.dy},{g.dx}) g{g.g_before}->{g.g_after}")
    for i, j in tr.recolours:
        parts.append(f"recolour c{A0.comps[i].color}->c{A1.comps[j].color} n{A0.comps[i].n}")
    rots = {(i, j) for i, j, _ in tr.inplace_rot}
    for i, j in tr.reshaped:
        tag = " (congruent, D4)" if (i, j) in rots else ""
        parts.append(f"reshape c{A0.comps[i].color} n{A0.comps[i].n}->n{A1.comps[j].n}{tag}")
    for i in tr.occluded:
        parts.append(f"occluded c{A0.comps[i].color} n{A0.comps[i].n}")
    for i in tr.annihilated:
        parts.append(f"annihilate c{A0.comps[i].color} n{A0.comps[i].n}")
    for j in tr.created:
        parts.append(f"create c{A1.comps[j].color} n{A1.comps[j].n}")
    return "; ".join(parts) if parts else "null"
