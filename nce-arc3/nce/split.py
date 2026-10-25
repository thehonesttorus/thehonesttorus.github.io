"""Downward refinement of strata under residuals (Section 2.5, two-projection graft) [D].

When a transition changes a strict subset of a component's cells, the component is re-filtered
by the Dirichlet form restricted to its changed and unchanged cells: the HN filtration of the
transition residual. The *unchanged* part is the stable stratum; the changed part that
transports exactly by a translation is a mover of its colour. A rule records the canonical
image of the unchanged remainder. At perception, a component of that colour that contains the
remainder's mask at exactly one position (exact sub-mask containment, found by
cross-correlation: second-order statistics, C11) is split into the remainder and the leftover
pieces, so later transitions see a plain gauge move of the leftover. Strata thereby earn their
status by stability under the observed correspondences.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .algebra import GRID, BaseAlgebra, canonical_form, d4_apply, d4_inverse, from_components, label4


def _mask_of(cells: np.ndarray) -> np.ndarray:
    m = np.zeros(GRID * GRID, dtype=bool)
    m[cells] = True
    return m.reshape(GRID, GRID)


def _containment_positions(big: np.ndarray, img: np.ndarray) -> list[tuple[int, int]]:
    """All bbox corners where `img` lies entirely inside `big`."""
    h, w = img.shape
    n = int(img.sum())
    if n == 0 or h > GRID or w > GRID:
        return []
    fa = np.fft.rfft2(img.astype(np.float32), s=big.shape)
    fb = np.fft.rfft2(big.astype(np.float32))
    corr = np.fft.irfft2(fb * np.conj(fa), s=big.shape)
    out = []
    for y, x in np.argwhere(corr > n - 0.5).tolist():
        if y + h <= GRID and x + w <= GRID and np.all(big[y:y + h, x:x + w][img]):
            out.append((y, x))
    return out


@dataclass
class SplitRules:
    remainders: dict[int, dict[bytes, np.ndarray]] = field(default_factory=dict)   # colour -> key -> canonical image of the stable part
    learned: int = 0

    def empty(self) -> bool:
        return not self.remainders

    # ---- learning from a transition --------------------------------------------------
    def learn(self, A0: BaseAlgebra, A1: BaseAlgebra) -> bool:
        """Learn stable-remainder types from reshaped pairs; True if a new rule was added."""
        added = False
        from .match import decompose
        tr = decompose(A0, A1)
        for i, j in tr.reshaped:
            a, b = A0.comps[i], A1.comps[j]
            if a.stratum == 0 and len(A0.strata) > 1:
                continue
            ma, mb = _mask_of(a.cells), _mask_of(b.cells)
            gone, new, stay = ma & ~mb, mb & ~ma, ma & mb
            if not gone.any() or not new.any() or not stay.any() or gone.sum() != new.sum():
                continue
            ys, xs = np.nonzero(gone); ys2, xs2 = np.nonzero(new)
            cg = gone[ys.min():ys.max() + 1, xs.min():xs.max() + 1]
            cn = new[ys2.min():ys2.max() + 1, xs2.min():xs2.max() + 1]
            if cg.shape != cn.shape or not np.array_equal(cg, cn):
                continue   # the changed part did not transport exactly
            sy, sx = np.nonzero(stay)
            crop = stay[sy.min():sy.max() + 1, sx.min():sx.max() + 1]
            key, g, canon = canonical_form(crop)
            table = self.remainders.setdefault(int(a.color), {})
            if key not in table:
                table[key] = canon
                added = True
                self.learned += 1
        return added

    # ---- application at perception --------------------------------------------------
    def apply(self, A: BaseAlgebra) -> BaseAlgebra:
        if not self.remainders:
            return A
        parts: list[tuple[int, np.ndarray]] = []
        changed = False
        for c in A.comps:
            table = self.remainders.get(c.color)
            if not table or (c.stratum == 0 and len(A.strata) > 1):
                parts.append(c); continue
            big = c.mask()
            done = False
            for key, canon in sorted(table.items(), key=lambda kv: -int(kv[1].sum())):
                n = int(canon.sum())
                if n >= c.n:
                    continue
                hits: list[np.ndarray] = []
                for g in range(8):
                    img = d4_apply(canon, d4_inverse(g))
                    for (y, x) in _containment_positions(big, img):
                        sub = np.zeros_like(big)
                        sub[y:y + img.shape[0], x:x + img.shape[1]] = img
                        if not any(np.array_equal(sub, h) for h in hits):
                            hits.append(sub)
                if len(hits) != 1:
                    continue           # ambiguous or absent: the rule does not apply here
                stay = hits[0]
                leftover = big & ~stay
                parts.append((c.color, np.flatnonzero(stay.ravel())))
                lab, k = label4(leftover)
                for t in range(1, k + 1):
                    parts.append((c.color, np.flatnonzero((lab == t).ravel())))
                changed = True
                done = True
                break
            if not done:
                parts.append(c)
        return from_components(A.frame, parts) if changed else A
