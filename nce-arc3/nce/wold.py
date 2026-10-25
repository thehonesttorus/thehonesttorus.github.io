"""Wold newborn decomposition on an increasing filtration (Theorem 1.4).

F_{<=l} = F_{<=l-1} (+) H_l, with H_l = ker D_l the newborn subspace and
Pi_new^(l) = Pi_{<=l} - Pi_{<=l-1}. The newborn component of a vector is what is orthogonal
to everything already explained; it is obtained by *signed subtraction* of the lower-grade
projection (Barrier I: positivity-preserving averaging cannot do this).

`WoldFiltration` maintains an orthonormal basis of the explained subspace incrementally
and returns, for any vector, its newborn part and norm. `residual` and
`transformed_residuals` implement the grade-0 version (observed minus predicted) and the
"make it linear before you match it" transformations of amendment item 2 / brief C11.
"""
from __future__ import annotations

from typing import Iterable, Optional

import numpy as np


class WoldFiltration:
    """Incremental orthonormal basis of the explained subspace, graded by insertion grade."""

    def __init__(self, dim: int, tol: float = 1e-9):
        self.dim = dim
        self.tol = tol
        self._basis = np.zeros((dim, 0))
        self._grades: list[int] = []

    @property
    def rank(self) -> int:
        return self._basis.shape[1]

    def project(self, x: np.ndarray) -> np.ndarray:
        """Pi_{<=l} x: the part already explained."""
        if self.rank == 0:
            return np.zeros_like(x, dtype=float)
        B = self._basis
        return B @ (B.T @ x)

    def newborn(self, x: np.ndarray) -> tuple[np.ndarray, float]:
        """(I - Pi) x and its norm: the genuinely new component."""
        r = np.asarray(x, dtype=float) - self.project(x)
        return r, float(np.linalg.norm(r))

    def absorb(self, x: np.ndarray, grade: int) -> bool:
        """Add the newborn direction of x to the explained subspace; True if it was new."""
        r, nr = self.newborn(x)
        if nr <= self.tol * max(1.0, float(np.linalg.norm(x))):
            return False
        r /= nr
        # re-orthogonalise once for numerical stability
        r -= self.project(r)
        r /= max(float(np.linalg.norm(r)), 1e-300)
        self._basis = np.concatenate([self._basis, r[:, None]], axis=1)
        self._grades.append(grade)
        return True

    def grade_projector(self, grade: int) -> np.ndarray:
        cols = [j for j, g in enumerate(self._grades) if g == grade]
        if not cols:
            return np.zeros((self.dim, self.dim))
        B = self._basis[:, cols]
        return B @ B.T

    def newborn_fraction(self, x: np.ndarray) -> float:
        """||(I - Pi) x||^2 / ||x||^2  (the newborn-ratio diagnostic of amendment item 11)."""
        n2 = float(np.dot(x, x))
        if n2 == 0.0:
            return 0.0
        _, nr = self.newborn(x)
        return nr * nr / n2


def residual(observed: np.ndarray, predicted: np.ndarray) -> np.ndarray:
    """Grade-0 newborn residual r = observed - predicted (signed subtraction)."""
    return np.asarray(observed, dtype=float) - np.asarray(predicted, dtype=float)


def transformed_residuals(r: np.ndarray, degree: int = 2) -> dict[str, np.ndarray]:
    """Label transformations that make degree-k structure linear ("make it linear before you match it").

    Returns the residual, its elementwise square, and (for degree >= 2) the pairwise
    products needed for a quadratic match. Amendment E7: a polynomial link has generative
    exponent <= 2, so matching on r and r^2 needs K ~ d samples instead of d^k.
    """
    r = np.asarray(r, dtype=float).ravel()
    out = {"r": r, "r2": r * r}
    if degree >= 2:
        out["rr"] = np.outer(r, r)[np.triu_indices(r.size)]
    return out


def hankel_rank_zero(residuals: Iterable[np.ndarray], atol: float = 0.0) -> bool:
    """Exact verification (T6): every stored residual is identically zero."""
    for r in residuals:
        if np.any(np.abs(np.asarray(r)) > atol):
            return False
    return True
