"""Kahler-Tyler moment map, Tyler fixed point, Forster-Barthe stability, Harder-Narasimhan flag.

Theory: Theorem 3.1 (i)-(iii) of the monograph (sound; Kent-Tyler 1988 for existence).

Given a measure nu on points s(x) in R^m \\ {0}, the potential
    Phi(G) = 1/2 log det G + m/2 * int log <s, G^{-1} s> dnu
is geodesically convex on positive-definite G; its moment map is
    mu_act(G) = I/m - int  s~ s~^T dnu,   s~ = G^{-1/2} s / |G^{-1/2} s|,
and a zero G* (Tyler's fixed point) exists and is unique up to scale iff the measure is
subspace-stable: nu(s in V) < dim V / m for every proper subspace V. When stability
fails, the gradient flow separates the eigenvalues of G exponentially across the
Harder-Narasimhan flag 0 = V_0 < V_1 < ... < V_r = R^m ordered by decreasing slope
mu_k = nu(s in V_k \\ V_{k-1}) / dim(V_k / V_{k-1}).

This module implements all three as plain numpy, with no assumptions about what the
points are. Scale is fixed by trace(G) = m (positive homothety is a gauge).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np


def _sym_inv_sqrt(G: np.ndarray) -> np.ndarray:
    w, U = np.linalg.eigh(G)
    w = np.clip(w, 1e-300, None)
    return (U / np.sqrt(w)) @ U.T


def moment_map(S: np.ndarray, G: np.ndarray, weights: Optional[np.ndarray] = None) -> np.ndarray:
    """mu_act(G) = I/m - sum_i w_i s~_i s~_i^T  (zero at the balanced metric)."""
    n, m = S.shape
    w = np.full(n, 1.0 / n) if weights is None else weights / weights.sum()
    T = S @ _sym_inv_sqrt(G)                      # rows G^{-1/2} s_i
    norms = np.linalg.norm(T, axis=1, keepdims=True)
    T = T / np.clip(norms, 1e-300, None)
    return np.eye(m) / m - (T * w[:, None]).T @ T


@dataclass
class TylerResult:
    G: np.ndarray                 # metric with trace m
    converged: bool
    residual: float               # Frobenius norm of the moment map at G
    iterations: int
    log_eigs: np.ndarray          # log eigenvalues of G (ascending), for HN diagnostics


def tyler_fixed_point(S: np.ndarray, weights: Optional[np.ndarray] = None, iters: int = 500,
                      tol: float = 1e-9, G0: Optional[np.ndarray] = None) -> TylerResult:
    """Tyler's fixed-point iteration  G <- m * sum_i w_i s_i s_i^T / <s_i, G^{-1} s_i>, trace-normalised.

    This is the alternating radial-normalise / whiten iteration (Corollary 5.2(i) of the
    monograph). Under subspace stability it converges to the unique balanced metric; under
    instability the eigenvalues of the trace-normalised iterate separate exponentially,
    which `hn_filtration` reads off.
    """
    n, m = S.shape
    w = np.full(n, 1.0 / n) if weights is None else weights / weights.sum()
    G = np.eye(m) if G0 is None else G0 / np.trace(G0) * m
    res = float("inf")
    for it in range(1, iters + 1):
        Ginv = np.linalg.pinv(G, hermitian=True)
        q = np.einsum("ij,jk,ik->i", S, Ginv, S)          # <s_i, G^{-1} s_i>
        q = np.clip(q, 1e-300, None)
        G_new = m * (S * (w / q)[:, None]).T @ S
        G_new = 0.5 * (G_new + G_new.T)
        tr = np.trace(G_new)
        if not np.isfinite(tr) or tr <= 0:
            break
        G_new *= m / tr
        res = float(np.linalg.norm(moment_map(S, G_new, w)))
        step = float(np.abs(G_new - G).max())
        G = G_new
        if step < tol:
            return TylerResult(G, True, res, it, np.log(np.clip(np.linalg.eigvalsh(G), 1e-300, None)))
    return TylerResult(G, False, res, iters, np.log(np.clip(np.linalg.eigvalsh(G), 1e-300, None)))


def whiten(S: np.ndarray, G: np.ndarray) -> np.ndarray:
    """Canonical radially-normalised states s~ = G^{-1/2} s / |G^{-1/2} s|  (points of S^{m-1})."""
    T = S @ _sym_inv_sqrt(G)
    return T / np.clip(np.linalg.norm(T, axis=1, keepdims=True), 1e-300, None)


@dataclass
class HNFlag:
    """Harder-Narasimhan flag read from the Tyler flow.

    blocks[k] is an orthonormal basis (columns) of V_k / V_{k-1}; slopes[k] is
    nu(s in V_k \\ V_{k-1}) / dim, decreasing in k. stable is True iff there is one block
    (the measure is subspace-stable and a balanced metric exists).
    """
    blocks: list[np.ndarray]
    slopes: list[float]
    masses: list[float]
    stable: bool


def hn_filtration(S: np.ndarray, weights: Optional[np.ndarray] = None, iters: int = 400,
                  gap: float = 2.0, atol: float = 1e-3) -> HNFlag:
    """Compute the HN flag by running the Tyler flow and clustering the log-eigenvalues of G.

    Under instability the eigenvalues of the trace-normalised iterate split into groups
    separated by gaps that grow linearly in log-scale; a gap larger than `gap` (nats)
    between consecutive sorted log-eigenvalues marks a block boundary. Points are
    assigned to the smallest V_k containing them (component outside V_k below atol),
    and slopes are mass / dimension. For a stable measure a single block is returned.
    """
    n, m = S.shape
    w = np.full(n, 1.0 / n) if weights is None else weights / weights.sum()
    res = tyler_fixed_point(S, w, iters=iters, tol=0.0)      # run the full flow
    evals, U = np.linalg.eigh(res.G)                            # ascending
    logs = np.log(np.clip(evals, 1e-300, None))
    # Descending order: largest eigenvalues of G are the directions the measure over-fills
    # (they are the top of the flag, V_1, with the largest slope).
    order = np.argsort(-evals)
    logs_d, U_d = logs[order], U[:, order]
    cuts = [0]
    for j in range(1, m):
        if logs_d[j - 1] - logs_d[j] > gap:
            cuts.append(j)
    cuts.append(m)
    blocks = [U_d[:, a:b] for a, b in zip(cuts, cuts[1:])]
    if len(blocks) == 1:
        return HNFlag(blocks, [1.0 / m], [1.0], True)
    # membership: smallest k with |P_{V_k}^perp s| < atol * |s|
    masses: list[float] = []
    assigned = np.zeros(n, dtype=bool)
    basis_so_far = np.zeros((m, 0))
    for B in blocks:
        basis_so_far = np.concatenate([basis_so_far, B], axis=1)
        proj = S @ basis_so_far @ basis_so_far.T
        outside = np.linalg.norm(S - proj, axis=1) / np.clip(np.linalg.norm(S, axis=1), 1e-300, None)
        inV = (outside < atol) & ~assigned
        masses.append(float(w[inV].sum()))
        assigned |= inV
    slopes = [mass / B.shape[1] for mass, B in zip(masses, blocks)]
    return HNFlag(blocks, slopes, masses, False)


def is_subspace_stable(S: np.ndarray, weights: Optional[np.ndarray] = None, iters: int = 300) -> bool:
    """Forster-Barthe / Kent-Tyler stability test via convergence of the Tyler flow."""
    r = tyler_fixed_point(S, weights, iters=iters)
    return r.converged and r.residual < 1e-6
