"""Bethe pair-belief estimator (prototype, numpy, float64).

State of a layer (pre-activations z): node beliefs (m, v, k3, k4 per neuron) and cross cumulants
C (pairs), S (full symmetric third-cumulant tensor; reference mode) or K = S[a,a,b] (factored modes).
The ReLU map expands every joint statistic of a = relu(z) around the product of the node beliefs:
E[prod f_a(z_a)] = E_prod[exp(D) prod f_a], D = sum over cross cumulants kappa_alpha d^alpha / alpha!,
and keeps the terms of order <= n^-1 (c ~ n^-1/2, K, T ~ n^-1) plus c^3, cK.
The linear map z' = a W is exact multilinear algebra.
"""
import numpy as np
from scipy.special import ndtr
from math import factorial

SQ2PI = np.sqrt(2 * np.pi)


def ladder(m, v, k3, k4, qmin=-8, qmax=4):
    """L_q = G_q + k3/6 G_{q-3} + k4/24 G_{q-4} for q in [qmin+6, qmax], G_q = E relu^q / q! (q>=0),
    G_{-k} = d^k/dm^k Phi.  Returns dict q -> array.  E_node[(relu^p)^{(r)}] = p! L_{p-r}."""
    s = np.sqrt(v)
    t = m / s
    phi = np.exp(-0.5 * t * t) / SQ2PI
    Phi = ndtr(t)
    G = {}
    M = [Phi, m * Phi + s * phi]
    for k in range(2, qmax + 1):
        M.append(m * M[k - 1] + (k - 1) * v * M[k - 2])
    for k in range(qmax + 1):
        G[k] = M[k] / factorial(k)
    # He_k(t) by recursion
    He = [np.ones_like(t), t]
    for k in range(2, -qmin + 2):
        He.append(t * He[k - 1] - (k - 1) * He[k - 2])
    for k in range(1, -qmin - 4 + 6):
        G[-k] = (-1) ** (k - 1) * He[k - 1] * phi / s ** k
    L = {}
    for q in range(qmin + 4, qmax + 1):
        L[q] = G[q] + k3 / 6 * G[q - 3] + k4 / 24 * G[q - 4]
    return L


def lam(L, p, r):
    return factorial(p) * L[p - r]


def node_moments(L):
    E1, E2, E3, E4 = L[1], 2 * L[2], 6 * L[3], 24 * L[4]
    mu = E1
    var = E2 - E1 ** 2
    k3 = E3 - 3 * E2 * E1 + 2 * E1 ** 3
    k4 = E4 - 4 * E3 * E1 - 3 * E2 ** 2 + 12 * E2 * E1 ** 2 - 6 * E1 ** 4
    return mu, var, k3, k4


def pair_delta(La, Lb, p, q, c, Kab, Kba, order3=True):
    """E[a^p b^q] - E a^p E b^q for a = relu(z_a), b = relu(z_b) (broadcast arrays: La index rows,
    Lb columns).  Terms: c, K, c^2, cK, c^3."""
    A = lambda r: lam(La, p, r)
    Bf = lambda r: lam(Lb, q, r)
    out = c * A(1) * Bf(1) + 0.5 * Kab * A(2) * Bf(1) + 0.5 * Kba * A(1) * Bf(2)
    out = out + 0.5 * c * c * A(2) * Bf(2)
    if order3:
        out = out + 0.5 * c * Kab * A(3) * Bf(2) + 0.5 * c * Kba * A(2) * Bf(3) + c ** 3 / 6 * A(3) * Bf(3)
    return out


def relu_map_full(m, C, S, k4, order3=True, keep_old=True):
    """Reference mode: full S (n,n,n).  Returns mu, Ca, Sa, k4a."""
    n = m.shape[0]
    v = np.diag(C).copy()
    k3 = np.einsum('aaa->a', S).copy()
    L = ladder(m, v, k3, k4)
    mu, var, k3a, k4a = node_moments(L)
    Lr = {q: x[:, None] for q, x in L.items()}
    Lc = {q: x[None, :] for q, x in L.items()}
    K = np.einsum('aab->ab', S)            # K_ab = S[a,a,b]
    c = C.copy(); np.fill_diagonal(c, 0.0)
    Kab = K.copy(); np.fill_diagonal(Kab, 0.0)
    Kba = Kab.T
    D11 = pair_delta(Lr, Lc, 1, 1, c, Kab, Kba, order3)
    D21 = pair_delta(Lr, Lc, 2, 1, c, Kab, Kba, order3)
    Ca = D11.copy(); np.fill_diagonal(Ca, var)
    Ka = D21 - 2 * mu[:, None] * D11; np.fill_diagonal(Ka, k3a)
    # distinct triples: T^z L0L0L0 + hubs
    L0, Lm1 = L[0], L[-1]
    Sa = np.zeros_like(S)
    if keep_old:
        Sa += S * (L0[:, None, None] * L0[None, :, None] * L0[None, None, :])
    # hub at b: c_ab c_bc L0a L-1b L0c, summed over the three centres
    u = c * L0[:, None]                    # u[a,b] = c_ab L0_a
    hub_b = np.einsum('ab,cb,b->abc', u, u, Lm1)  # centre b
    Sa += hub_b + hub_b.transpose(1, 0, 2) + hub_b.transpose(0, 2, 1)  # centre b, a, c
    # overwrite the (a,a,b) patterns and diagonal with the pair/node values
    idx = np.arange(n)
    Sa[idx, idx, :] = Ka
    Sa[idx, :, idx] = Ka
    Sa[:, idx, idx] = Ka.T
    Sa[idx, idx, idx] = k3a
    return mu, Ca, Sa, k4a


def linear_full(mu, Ca, Sa, k4a, W):
    m = mu @ W
    C = W.T @ Ca @ W
    S = np.tensordot(Sa, W, axes=([2], [0]))
    S = np.tensordot(S, W, axes=([1], [0])).transpose(0, 2, 1)
    S = np.tensordot(W, S, axes=([0], [0]))
    k4 = (W ** 4).T @ k4a
    return m, C, S, k4


def estimate_full(Ws, order3=True, keep_old=True):
    """Ws: (L, n, n) in the x @ W convention.  Returns (L, n) means."""
    Ls, n, _ = Ws.shape
    W = Ws[0].astype(np.float64)
    m = np.zeros(n); C = W.T @ W; S = np.zeros((n, n, n)); k4 = np.zeros(n)
    out = []
    for l in range(Ls):
        mu, Ca, Sa, k4a = relu_map_full(m, C, S, k4, order3, keep_old)
        out.append(mu)
        if l + 1 < Ls:
            m, C, S, k4 = linear_full(mu, Ca, Sa, k4a, Ws[l + 1].astype(np.float64))
    return np.array(out)


def estimate_gauss(Ws):
    """Baseline: Gaussian closure with exact bivariate ReLU kernel (arc-cosine with means)."""
    from scipy.stats import multivariate_normal  # noqa (not used; closed form below via c-series)
    raise NotImplementedError


def estimate_gauss(Ws, R=12):
    """Baseline: Gaussian (covariance) closure; bivariate ReLU kernel by the Mehler series to order R."""
    Ls, n, _ = Ws.shape
    W = Ws[0].astype(np.float64)
    m = np.zeros(n); C = W.T @ W
    out = []
    for l in range(Ls):
        v = np.diag(C).copy()
        L = ladder(m, v, 0 * v, 0 * v, qmin=-R - 6, qmax=4)
        mu, var, _, _ = node_moments(L)
        c = C.copy(); np.fill_diagonal(c, 0)
        Ca = np.zeros_like(C); cr = np.ones_like(C); f = 1.0
        for r in range(1, R + 1):
            cr = cr * c; f *= r
            Ca += cr / f * np.outer(L[1 - r], L[1 - r])
        np.fill_diagonal(Ca, var)
        out.append(mu)
        if l + 1 < Ls:
            W = Ws[l + 1].astype(np.float64)
            m = mu @ W; C = W.T @ Ca @ W
    return np.array(out)


def mc_truth(Ws, N, batch=1_000_000, seed=0):
    rng = np.random.default_rng(seed)
    Ls, n, _ = Ws.shape
    W32 = Ws.astype(np.float32)
    acc = np.zeros((Ls, n)); acc2 = np.zeros((Ls, n))
    done = 0
    while done < N:
        b = min(batch, N - done)
        a = rng.standard_normal((b, n), dtype=np.float32)
        for l in range(Ls):
            a = np.maximum(a @ W32[l], 0)
            acc[l] += a.sum(0, dtype=np.float64); acc2[l] += (a.astype(np.float64) ** 2).sum(0)
        done += b
    mean = acc / N
    var = acc2 / N - mean ** 2
    return mean, var / N   # means and their MC variance
