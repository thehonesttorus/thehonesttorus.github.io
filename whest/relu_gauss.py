"""Exact Gaussian moment maps of ReLU, and the Gaussian (NNGP-type, quenched-weight) closure chain.

Conventions. Weights W[l] have shape (n_out, n_in) in column convention: z_l = W[l] @ h_{l-1}, h_l = relu(z_l),
h_0 = x ~ N(0, I). The truth file stores the post-activation means m[l-1] = E[h_l] for l = 1..L.

Hermite coefficients of r_a(u) = relu(u + a) for u ~ N(0,1), with He_k the probabilists' Hermite polynomials:
    d_0 = phi(a) + a Phi(a),  d_1 = Phi(a),  d_k = phi(a) He_{k-2}(-a)  (k >= 2),
so that for standardised (u, v) with correlation rho,  E[r_a(u) r_b(v)] = sum_k rho^k d_k(a) d_k(b) / k!.
"""
import numpy as np
from scipy.special import ndtr

SQ2PI = np.sqrt(2 * np.pi)


def phi(a):
    return np.exp(-0.5 * a * a) / SQ2PI


def hermite_relu(alpha, K):
    """d_k(alpha) for k = 0..K, shape (K+1, *alpha.shape)."""
    alpha = np.asarray(alpha, dtype=np.float64)
    P = ndtr(alpha); p = phi(alpha)
    d = np.empty((K + 1,) + alpha.shape)
    d[0] = p + alpha * P
    if K >= 1:
        d[1] = P
    x = -alpha; Hm, Hc = np.zeros_like(x), np.ones_like(x)   # He_{j-1}, He_j with j = 0
    for k in range(2, K + 1):
        j = k - 2
        if j > 0:
            Hm, Hc = Hc, x * Hc - (j - 1) * Hm if j > 1 else x * Hc  # He_j = x He_{j-1} - (j-1) He_{j-2}
        d[k] = p * Hc
    return d


def relu_moments(mu, C, K=10):
    """Mean and covariance of h = relu(z), z ~ N(mu, C), by the Hermite series in the off-diagonal correlations and
    the exact diagonal. Returns (m, Kh, sigma, alpha)."""
    var = np.clip(np.diag(C), 1e-30, None); sigma = np.sqrt(var); alpha = mu / sigma
    d = hermite_relu(alpha, K)
    m = sigma * d[0]
    rho = C / np.outer(sigma, sigma)
    np.fill_diagonal(rho, 0.0)
    Kh = np.zeros_like(C); rk = np.ones_like(rho); fact = 1.0
    for k in range(1, K + 1):
        rk = rk * rho; fact *= k
        Kh += np.outer(d[k], d[k]) * rk / fact
    Kh *= np.outer(sigma, sigma)
    second = var * ((1 + alpha * alpha) * ndtr(alpha) + alpha * phi(alpha))   # E[relu(z)^2]
    np.fill_diagonal(Kh, second - m * m)
    return m, Kh, sigma, alpha


def gauss_chain(W, K=10, record=None):
    """Gaussian closure: propagate (mu, C) of the pre-activations layer by layer. Returns the post-activation means
    of every layer, shape (L, n). `record`, if a dict, receives per-layer mu, C, sigma, alpha, m, Kh."""
    L, n, n_in = W.shape
    mu = np.zeros(n_in); C = np.eye(n_in); out = np.empty((L, n))
    for l in range(L):
        Wl = W[l].astype(np.float64)
        mu_z = Wl @ mu; C_z = Wl @ C @ Wl.T
        m, Kh, sigma, alpha = relu_moments(mu_z, C_z, K)
        out[l] = m
        if record is not None:
            record[l] = dict(mu=mu_z, C=C_z, sigma=sigma, alpha=alpha, m=m, Kh=Kh)
        mu, C = m, Kh
    return out
