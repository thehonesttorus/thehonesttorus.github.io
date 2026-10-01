"""Gaussian covariance closures at n = 1024 (numpy, float64), with hooks for linear-response probes.

State per layer: pre-activation mean mu_l, covariance S_l = C(z_l); post-activation mean m_l and covariance C_l = C(a_l).
Closures for C_l from (mu_l, S_l):
  'lin'   : off-diagonal Phi_a Phi_b S_ab (the bench baseline, linearised cross-covariance)
  'exact' : the exact bivariate-Gaussian ReLU covariance,
            Cov(relu z_a, relu z_b) = c Phi_a Phi_b + int_0^c (c - t) p_t(0, 0) dt,   c = S_ab,
            p_t(0,0) = bivariate normal density of (z_a, z_b) at the origin with covariance t (d^2/dc^2 E[f g] = E[f'' g''])
            evaluated by Gauss-Legendre in t (exact for the Gaussian law up to quadrature error).
"""
import numpy as np
from scipy.special import ndtr

SQ2PI = np.sqrt(2 * np.pi)


def marg(mu, v):
    v = np.maximum(v, 1e-300); s = np.sqrt(v); al = mu / s
    P = ndtr(al); p = np.exp(-0.5 * al * al) / SQ2PI
    m = mu * P + s * p
    sec = (mu * mu + v) * P + mu * s * p
    return m, sec - m * m, P, p, s


_GL = {}


def gl(k):
    if k not in _GL:
        x, w = np.polynomial.legendre.leggauss(k)
        _GL[k] = (0.5 * (x + 1), 0.5 * w)
    return _GL[k]


def cov_exact(mu, S, P, k=24):
    """Off-diagonal exact Gaussian ReLU covariance (diagonal filled by caller)."""
    v = np.diag(S).copy()
    out = S * P[:, None] * P[None, :]
    nodes, wts = gl(k)
    va = v[:, None]; vb = v[None, :]
    mua = mu[:, None]; mub = mu[None, :]
    vv = va * vb
    num0 = mua * mua * vb + mub * mub * va
    mm = mua * mub
    acc = np.zeros_like(S)
    for s_, w_ in zip(nodes, wts):
        t = S * s_
        det = np.maximum(vv - t * t, 1e-12 * vv)
        Q = (num0 - 2 * t * mm) / det
        acc += w_ * (1 - s_) * np.exp(-0.5 * Q) / (2 * np.pi * np.sqrt(det))
    out += S * S * acc
    return out


def closure_C(mu, S, kind="exact"):
    m, var, P, p, s = marg(mu, np.diag(S))
    if kind == "lin":
        C = S * P[:, None] * P[None, :]
    else:
        C = cov_exact(mu, S, P)
    np.fill_diagonal(C, var)
    return m, C, P, p, s


def run(W, kind="exact", hook=None, start=None):
    """Forward Gaussian closure. hook(l, m, C) -> (m, C) may perturb the post-activation state of layer l.
    start = (l0, m, C): resume from the post-activation state of layer l0 - 1 (for probes). Returns (L, n) means and
    the list of post-activation (m, C)."""
    W = W.astype(np.float64) if W.dtype != np.float64 else W
    L, n, _ = W.shape
    outs, states = [], []
    if start is None:
        l0 = 0; m = None; C = None
    else:
        l0, m, C = start
    for l in range(l0, L):
        if l == 0:
            mu = np.zeros(n); S = W[0].T @ W[0]
        else:
            mu = m @ W[l]; S = W[l].T @ C @ W[l]
        m, C, P, p, s = closure_C(mu, S, kind)
        if hook is not None:
            m, C = hook(l, m, C)
        outs.append(m); states.append((m, C))
    return np.stack(outs), states


def predict_exact(W):
    return run(W, "exact")[0]


def predict_lin(W):
    return run(W, "lin")[0]
