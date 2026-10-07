# Exact checks of the mean-preserving radial tangent (note XXXII section 4).   python check_radial.py
# Y = S A with S = 1 + eps independent of A, E eps = 0, E eps^2 = v, E eps^3 = 0. To first order in v:
#   d mu = 0,  d Sigma = v (Sigma + mu mu^T),  d kappa3 = v (3 kappa3 + 2 Sym3(mu, Sigma)),
#   d kappa4 = v (6 kappa4 + 4 Sym3(Sigma, Sigma) + 3 Sym4(mu, kappa3)).
# (a) checked on an arbitrary finite law of A in R^2 by exact enumeration with a symmetric two-point eps;
# (b) the Gram-Charlier mean map E relu = mu Phi + sigma phi - kappa3 a phi/(6 sigma^2) + kappa4 (a^2-1) phi/(24 sigma^3)
#     has zero first-order response to this tangent at a Gaussian point (the cancellation is exact in alpha).
import itertools, numpy as np
rng = np.random.default_rng(3)
pts = rng.standard_normal((7, 2)) + np.array([0.6, -0.3]); w = rng.random(7); w /= w.sum()   # finite law of A


def cumulants(P, wt):
    mu = wt @ P; X = P - mu
    S2 = np.einsum("k,ki,kj->ij", wt, X, X)
    C3 = np.einsum("k,ki,kj,kl->ijl", wt, X, X, X)
    M4 = np.einsum("k,ki,kj,kl,km->ijlm", wt, X, X, X, X)
    K4 = M4 - (np.einsum("ij,lm->ijlm", S2, S2) + np.einsum("il,jm->ijlm", S2, S2) + np.einsum("im,jl->ijlm", S2, S2))
    return mu, S2, C3, K4


mu, Sg, C3, K4 = cumulants(pts, w)
def scaled(v):
    e = np.sqrt(v)
    return cumulants(np.vstack([(1 + e) * pts, (1 - e) * pts]), np.concatenate([w, w]) / 2)
h = 1e-4   # one-sided difference: exact for the linear terms, an O(v) residual from E eps^4 = v^2 in kappa4 (halves with h)
num = [(a - b) / h for a, b in zip(scaled(h), scaled(0.0))]
Sym3_mS = np.einsum("i,jk->ijk", mu, Sg) + np.einsum("j,ik->ijk", mu, Sg) + np.einsum("k,ij->ijk", mu, Sg)
Sym3_SS = np.einsum("ij,kl->ijkl", Sg, Sg) + np.einsum("ik,jl->ijkl", Sg, Sg) + np.einsum("il,jk->ijkl", Sg, Sg)
Sym4_mC = (np.einsum("i,jkl->ijkl", mu, C3) + np.einsum("j,ikl->ijkl", mu, C3) + np.einsum("k,ijl->ijkl", mu, C3)
           + np.einsum("l,ijk->ijkl", mu, C3))
pred = [np.zeros(2), Sg + np.outer(mu, mu), 3 * C3 + 2 * Sym3_mS, 6 * K4 + 4 * Sym3_SS + 3 * Sym4_mC]
for nm, a, b in zip(("d mu", "d Sigma", "d kappa3", "d kappa4"), num, pred):
    print(f"(a) {nm}: max |finite difference - formula| = {np.abs(a - b).max():.2e} (scale {np.abs(b).max():.2f}; O(v) step)")
# (b) Gram-Charlier first-order response along the tangent at a Gaussian point (kappa3 = kappa4 = 0)
from math import erf, exp, pi, sqrt
for (m, s) in [(0.0, 1.0), (0.7, 1.3), (-1.1, 0.6), (2.0, 0.9)]:
    a = m / s; ph = exp(-a * a / 2) / sqrt(2 * pi)
    dvar = s * s + m * m; dk3 = 2 * 3 * m * s * s; dk4 = 4 * 3 * s ** 4          # per unit v
    dmean = ph * dvar / (2 * s) - dk3 * a * ph / (6 * s * s) + dk4 * (a * a - 1) * ph / (24 * s ** 3)
    print(f"(b) mu {m:+.1f} sigma {s:.1f}: first-order mean response along the radial tangent = {dmean:+.2e}")
