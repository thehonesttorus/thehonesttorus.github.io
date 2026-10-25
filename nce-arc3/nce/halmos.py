"""Two-projection (Halmos) geometry and Kato transport utilities.

Theory: Theorem 4.1(i) (minimum-norm interpolation and the tangent tax; sound) and
Theorem 4.1(ii) with amendment E10 (||[P,Q]||_HS^2 = 1/2 Tr sin^2 2Theta = 2 Tr(P Q P^perp Q);
the Kato speed limit |dA/dt| <= ||G||_HS ||[P,Q]||_HS is attained for G proportional to [P,Q]).

Subspaces are given by matrices whose columns span them (any basis); orthonormal bases
are computed internally.
"""
from __future__ import annotations

import numpy as np


def orthonormal_basis(A: np.ndarray, tol: float = 1e-10) -> np.ndarray:
    """Orthonormal basis of the column span of A (empty (n,0) matrix for A = 0)."""
    if A.ndim == 1:
        A = A[:, None]
    if A.shape[1] == 0:
        return A
    U, s, _ = np.linalg.svd(A, full_matrices=False)
    r = int((s > tol * max(1.0, s[0] if s.size else 1.0)).sum())
    return U[:, :r]


def projector(A: np.ndarray) -> np.ndarray:
    B = orthonormal_basis(A)
    return B @ B.T


def principal_angles(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Principal angles (radians, ascending) between span(A) and span(B)."""
    QA, QB = orthonormal_basis(A), orthonormal_basis(B)
    if QA.shape[1] == 0 or QB.shape[1] == 0:
        return np.zeros(0)
    s = np.linalg.svd(QA.T @ QB, compute_uv=False)
    return np.arccos(np.clip(s, -1.0, 1.0))


def alignment(P: np.ndarray, Q: np.ndarray) -> float:
    """A = Tr(P Q): the subspace alignment functional of Theorem 4.1(ii)."""
    return float(np.trace(P @ Q))


def commutator_hs(P: np.ndarray, Q: np.ndarray) -> float:
    """||[P, Q]||_HS. Equals sqrt(2 Tr(P Q P^perp Q)) = sqrt(1/2 Tr sin^2 2Theta) (amendment E10)."""
    C = P @ Q - Q @ P
    return float(np.sqrt(np.sum(C * C)))


def kato_generator(P: np.ndarray, Pdot: np.ndarray) -> np.ndarray:
    """G = [Pdot, P]; satisfies Pdot = [G, P] and is off-diagonal w.r.t. P (Theorem 4.1(ii))."""
    return Pdot @ P - P @ Pdot


def kato_speed_limit(P: np.ndarray, Q: np.ndarray, G: np.ndarray) -> tuple[float, float]:
    """(|dA/dt|, bound) with dA/dt = Tr(G [P,Q]) and bound = ||G||_HS ||[P,Q]||_HS."""
    C = P @ Q - Q @ P
    rate = abs(float(np.trace(G @ C)))
    return rate, float(np.sqrt(np.sum(G * G)) * np.sqrt(np.sum(C * C)))


def min_norm_interpolant(P_basis: np.ndarray, Q_basis: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, float]:
    """Unique minimum-norm f in span(P_basis) with Q f = y, for y in span(Q_basis).

    Returns (f, tangent_tax) where ||f||^2 = ||y||^2 + tangent_tax and
    tangent_tax = ||tan Theta(P,Q) y||^2 >= 0 (Theorem 4.1(i)). Requires Q P H = Q H on the
    relevant part (all principal angles < pi/2 along y); otherwise the least-squares
    solution is returned and the identity need not hold.
    """
    P = projector(P_basis)
    Q = projector(Q_basis)
    QPQ = Q @ P @ Q
    # f = P Q (Q P Q|_{QH})^{-1} y, computed by pseudo-inverse on QH
    f = P @ Q @ np.linalg.pinv(QPQ, rcond=1e-12) @ y
    tax = float(f @ f - y @ y)
    return f, tax


def tangent_tax(P_basis: np.ndarray, y: np.ndarray) -> float:
    """||tan Theta(P, span y) y||^2 for a single target vector y (Q = projector onto y)."""
    _, tax = min_norm_interpolant(P_basis, y[:, None], y)
    return tax


def disagreement(predictions: np.ndarray) -> float:
    """Commutator-type disagreement among K predicted vectors (rows).

    For rank-one projections P_i onto the predictions, sum_{i<j} ||[P_i, P_j]||_HS^2 =
    sum_{i<j} 2 cos^2 theta_ij sin^2 theta_ij; this is the Kato criterion's quantity for a
    coherent set of hypotheses (T3). Zero vectors are ignored.
    """
    V = np.asarray(predictions, dtype=float)
    norms = np.linalg.norm(V, axis=1)
    V = V[norms > 0] / norms[norms > 0, None]
    if V.shape[0] < 2:
        return 0.0
    C = np.clip(V @ V.T, -1.0, 1.0)
    c2 = C * C
    iu = np.triu_indices(V.shape[0], k=1)
    return float(np.sum(2.0 * c2[iu] * (1.0 - c2[iu])))


def spread(predictions: np.ndarray) -> float:
    """Total variance of a set of predicted vectors (rows): sum_i ||p_i - mean||^2 / K.

    A positive-definite alternative to `disagreement` that does not vanish when
    predictions are antipodal; used when predictions are not naturally rays.
    """
    V = np.asarray(predictions, dtype=float)
    if V.shape[0] < 2:
        return 0.0
    return float(np.mean(np.sum((V - V.mean(axis=0)) ** 2, axis=1)))
