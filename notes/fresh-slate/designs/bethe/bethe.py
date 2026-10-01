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


# ------------------------------------------------------------------ production (factored) mode ---
def contract(k3, Kt, c, L0, Lm1, A, B):
    """S[A_j, A_j, B_k] for S = diag(k3) + (2,1)-patterns(Kt, zero diagonal) + hubs(c, L0, Lm1).
    Five n^3 products when A is B."""
    KA = Kt @ A
    KB = KA if B is A else Kt @ B
    X = k3[:, None] * A * A + 2 * A * KA
    out = None
    if c is not None:
        UA = c @ (L0[:, None] * A)
        UB = UA if B is A else c @ (L0[:, None] * B)
        X = X + Lm1[:, None] * UA * UA
        out = 2 * (Lm1[:, None] * A * UA).T @ UB
    K = X.T @ B + (A * A).T @ KB
    return K if out is None else K + out


def relu_map_pairs(m, C, K, k4, order3=True):
    v = np.diag(C).copy(); k3 = np.diag(K).copy()
    L = ladder(m, v, k3, k4)
    mu, var, k3a, k4a = node_moments(L)
    Lr = {q: x[:, None] for q, x in L.items()}
    Lc = {q: x[None, :] for q, x in L.items()}
    c = C.copy(); np.fill_diagonal(c, 0.0)
    Kab = K.copy(); np.fill_diagonal(Kab, 0.0)
    D11 = pair_delta(Lr, Lc, 1, 1, c, Kab, Kab.T, order3)
    D21 = pair_delta(Lr, Lc, 2, 1, c, Kab, Kab.T, order3)
    Ca = D11.copy(); np.fill_diagonal(Ca, var)
    Ka = D21 - 2 * mu[:, None] * D11; np.fill_diagonal(Ka, 0.0)
    return mu, Ca, Ka, k3a, k4a, c, L[0], L[-1]


def estimate_fact(Ws, old=1, hubs=True, order3=True):
    """Production mode: state (m, C, K, k4); distinct triples generated by hubs; old triple content
    carried to age `old` (0 or 1) through the gated two-step propagator."""
    Ls, n, _ = Ws.shape
    W = Ws[0].astype(np.float64)
    m = np.zeros(n); C = W.T @ W; K = np.zeros((n, n)); k4 = np.zeros(n)
    out = []; prev = None
    for l in range(Ls):
        mu, Ca, Ka, k3a, k4a, c, L0, Lm1 = relu_map_pairs(m, C, K, k4, order3)
        out.append(mu)
        if l + 1 == Ls:
            break
        Wn = Ws[l + 1].astype(np.float64)
        M = c * c * (L0 ** 2)[:, None] * Lm1[None, :] if hubs else 0.0   # hub value on (a,a,b)
        Kt = Ka - M
        spec = (k3a, Kt, c if hubs else None, L0, Lm1)
        Kn = contract(*spec, Wn, Wn)
        if old and prev is not None:
            pspec, Wl = prev
            P = Wl @ (L0[:, None] * Wn)                     # gated two-step propagator
            Kz = K.copy(); k3z = np.diag(K).copy(); np.fill_diagonal(Kz, 0.0)
            full = contract(*pspec, P, P)
            pat = contract(k3z * L0 ** 3, Kz * (L0 ** 2)[:, None] * L0[None, :], None, L0, Lm1, Wn, Wn)
            Kn = Kn + full - pat
        prev = (spec, Wn)
        m = mu @ Wn
        C = Wn.T @ Ca @ Wn
        K = Kn
        k4 = (Wn ** 4).T @ k4a
    return np.array(out)


def estimate_tree(Ws):
    """Pure node-belief Bethe (lift limit): independent inputs, cumulants by W^{o r}."""
    Ls, n, _ = Ws.shape
    W = Ws[0].astype(np.float64)
    m = np.zeros(n); v = (W ** 2).sum(0); k3 = np.zeros(n); k4 = np.zeros(n)
    out = []
    for l in range(Ls):
        L = ladder(m, v, k3, k4)
        mu, var, k3a, k4a = node_moments(L)
        out.append(mu)
        if l + 1 == Ls:
            break
        W = Ws[l + 1].astype(np.float64)
        m = mu @ W; v = var @ W ** 2; k3 = k3a @ W ** 3; k4 = k4a @ W ** 4
    return np.array(out)


# ------------------------------------------------- v1: exact bivariate-Gaussian edge beliefs -------
_GLX, _GLW = np.polynomial.legendre.leggauss(40)


def bvn_orthant(ta, tb, rho):
    """P(z_a > 0, z_b > 0) = Phi2(ta, tb; rho) by Gauss-Legendre in theta, rho = sin(theta)."""
    th1 = np.arcsin(rho)
    out = ndtr(ta) * ndtr(tb)
    acc = np.zeros(np.broadcast(ta, tb, rho).shape)
    hk = ta * tb; hh = 0.5 * (ta * ta + tb * tb)
    for x, w in zip(_GLX, _GLW):
        th = 0.5 * th1 * (x + 1)
        s = np.sin(th); c2 = np.cos(th) ** 2
        acc += w * np.exp(-(hh - hk * s) / c2)
    return out + acc * 0.5 * th1 / (2 * np.pi)


def uni_G(m, v, qmin, qmax):
    """Gaussian ladder G_q, q in [qmin, qmax]: G_q = E relu^q / q! (q >= 0), G_{-k} = d^k Phi / dm^k."""
    s = np.sqrt(v); t = m / s
    phi = np.exp(-0.5 * t * t) / SQ2PI; Phi = ndtr(t)
    G = {}
    M = [Phi, m * Phi + s * phi]
    for k in range(2, max(qmax, 1) + 1):
        M.append(m * M[k - 1] + (k - 1) * v * M[k - 2])
    for k in range(0, qmax + 1):
        G[k] = M[k] / factorial(k)
    He = [np.ones_like(t), t]
    for k in range(2, -qmin + 1):
        He.append(t * He[k - 1] - (k - 1) * He[k - 2])
    for k in range(1, -qmin + 1):
        G[-k] = (-1) ** (k - 1) * He[k - 1] * phi / s ** k
    return G


def edge_table(m, v, c, pmax=2, kmin=-4):
    """G_{ij} = E_G2[relu_a^i relu_b^j]/(i! j!) with negative indices = mean derivatives, for the
    bivariate Gaussian of every pair (a row, b column).  i, j in [kmin, pmax]."""
    n = m.shape[0]
    ma, mb = m[:, None], m[None, :]
    va, vb = v[:, None], v[None, :]
    rho = np.clip(c / np.sqrt(va * vb), -0.999999, 0.999999)
    c = rho * np.sqrt(va * vb)
    Ga = uni_G(ma, va, kmin - 1 + kmin, pmax)
    Gb = uni_G(mb, vb, kmin - 1 + kmin, pmax)
    vba = np.maximum(vb - c * c / va, 1e-12 * vb); mba = mb - c * ma / va
    vab = np.maximum(va - c * c / vb, 1e-12 * va); mab = ma - c * mb / vb
    gba = uni_G(mba, vba, kmin + kmin, pmax)
    gab = uni_G(mab, vab, kmin + kmin, pmax)
    from math import comb
    T = {}
    for k in range(1, -kmin + 1):
        for j in range(kmin, pmax + 1):
            T[(-k, j)] = sum(comb(k - 1, i) * Ga[-k + i] * (-c / va) ** i * gba[j - i] for i in range(k))
            if j >= 0:
                T[(j, -k)] = sum(comb(k - 1, i) * Gb[-k + i] * (-c / vb) ** i * gab[j - i] for i in range(k))
    T[(0, 0)] = bvn_orthant(ma / np.sqrt(va), mb / np.sqrt(vb), rho)
    for q in range(1, pmax + 1):
        T[(0, q)] = (mb * T[(0, q - 1)] + vb * T[(0, q - 2)] + c * T[(-1, q - 1)]) / q if q >= 2 else \
            (mb * T[(0, 0)] + vb * T[(0, -1)] + c * T[(-1, 0)])
    for p in range(1, pmax + 1):
        for q in range(0, pmax + 1):
            T[(p, q)] = (ma * T[(p - 1, q)] + va * T[(p - 2, q)] + c * T[(p - 1, q - 1)]) / p
    return T


def pair_moment(T, p, q, k3a, k3b, k4a, k4b, Kab, Kba, Q=None, Rab=None, Rba=None):
    g = lambda i, j: T[(i, j)]
    out = (g(p, q) + k3a / 6 * g(p - 3, q) + k3b / 6 * g(p, q - 3) + k4a / 24 * g(p - 4, q)
           + k4b / 24 * g(p, q - 4) + Kab / 2 * g(p - 2, q - 1) + Kba / 2 * g(p - 1, q - 2))
    if Q is not None:
        out = out + Q / 4 * g(p - 2, q - 2) + Rab / 6 * g(p - 3, q - 1) + Rba / 6 * g(p - 1, q - 3)
    return factorial(p) * factorial(q) * out


def relu_map_edges(m, C, K, k4):
    v = np.diag(C).copy(); k3 = np.diag(K).copy()
    L = ladder(m, v, k3, k4)
    mu, var, k3a, k4a = node_moments(L)
    c = C.copy(); np.fill_diagonal(c, 0.0)
    Kab = K.copy(); np.fill_diagonal(Kab, 0.0)
    T = edge_table(m, v, c)
    r, cc = (lambda x: x[:, None]), (lambda x: x[None, :])
    args = (r(k3), cc(k3), r(k4), cc(k4), Kab, Kab.T)
    E11 = pair_moment(T, 1, 1, *args)
    E21 = pair_moment(T, 2, 1, *args)
    E2 = 2 * L[2]
    Ca = E11 - np.outer(mu, mu)
    Ka = E21 - np.outer(E2, mu) - 2 * mu[:, None] * Ca
    np.fill_diagonal(Ca, var); np.fill_diagonal(Ka, 0.0)
    return mu, Ca, Ka, k3a, k4a, c, L[0], L[-1]


def estimate_edge(Ws, old=1, hubs=True):
    """v1: as estimate_fact but with exact bivariate-Gaussian edge beliefs (all orders in c)."""
    Ls, n, _ = Ws.shape
    W = Ws[0].astype(np.float64)
    m = np.zeros(n); C = W.T @ W; K = np.zeros((n, n)); k4 = np.zeros(n)
    out = []; prev = None
    for l in range(Ls):
        mu, Ca, Ka, k3a, k4a, c, L0, Lm1 = relu_map_edges(m, C, K, k4)
        out.append(mu)
        if l + 1 == Ls:
            break
        Wn = Ws[l + 1].astype(np.float64)
        M = c * c * (L0 ** 2)[:, None] * Lm1[None, :] if hubs else 0.0
        spec = (k3a, Ka - M, c if hubs else None, L0, Lm1)
        Kn = contract(*spec, Wn, Wn)
        if old and prev is not None:
            pspec, Wl = prev
            P = Wl @ (L0[:, None] * Wn)
            Kz = K.copy(); k3z = np.diag(K).copy(); np.fill_diagonal(Kz, 0.0)
            Kn = Kn + contract(*pspec, P, P) - contract(k3z * L0 ** 3, Kz * (L0 ** 2)[:, None] * L0[None, :],
                                                      None, L0, Lm1, Wn, Wn)
        prev = (spec, Wn)
        m = mu @ Wn; C = Wn.T @ Ca @ Wn; K = Kn; k4 = (Wn ** 4).T @ k4a
    return np.array(out)


# ------------------------------------------------- v2: node kappa4 with two-site and tree terms ----
def pair_cumulants4(T, mu, E2, E3, args):
    """Q_ab = kappa(a,a,b,b), R_ab = kappa(a,a,a,b) of a = relu(z) from the edge table."""
    E11 = pair_moment(T, 1, 1, *args); E21 = pair_moment(T, 2, 1, *args)
    E22 = pair_moment(T, 2, 2, *args); E31 = pair_moment(T, 3, 1, *args)
    E12 = E21.T
    ma, mb = mu[:, None], mu[None, :]
    Ea = {0: 1.0, 1: ma, 2: E2[:, None], 3: E3[:, None]}
    Eb = {0: 1.0, 1: mb, 2: E2[None, :]}
    Eab = {(1, 1): E11, (2, 1): E21, (1, 2): E12, (2, 2): E22, (3, 1): E31}
    def raw(i, j):
        if i == 0: return Eb[j]
        if j == 0: return Ea[i]
        return Eab[(i, j)]
    from math import comb
    def central(p, q):
        return sum(comb(p, i) * comb(q, j) * (-ma) ** (p - i) * (-mb) ** (q - j) * raw(i, j)
                   for i in range(p + 1) for j in range(q + 1))
    c11 = central(1, 1); c20 = central(2, 0); c02 = central(0, 2)
    R = central(3, 1) - 3 * c20 * c11
    Q = central(2, 2) - c20 * c02 - 2 * c11 ** 2
    np.fill_diagonal(R, 0.0); np.fill_diagonal(Q, 0.0)
    return Q, R


def k4_next(W, k4a, Q, R, c, L, mu):
    """kappa4 of z_j = sum_a W_aj a_a: node + (3,1) + (2,2) patterns + generated (2,1,1) and
    (1,1,1,1) tree terms (with the leading coincidence subtractions)."""
    L0, Lm1, Lm2 = L[0], L[-1], L[-2]
    W2 = W * W
    out = (W2 * W2).T @ k4a
    out += 4 * np.einsum('aj,aj->j', W2 * W, R @ W)
    out += 3 * np.einsum('aj,aj->j', W2, Q @ W2)
    X = L0[:, None] * W                     # x_ab = c_ab L0_b W_b
    U = c @ X
    c2 = c * c
    s = c2 @ (X * X)                        # sum_b x_ab^2
    h4 = 2 * (L0 * (1 - L0) - mu * Lm1)     # kappa(a,a,.,.) hub coefficient
    out += 6 * np.einsum('aj,aj->j', W2 * h4[:, None], U * U - s)
    kk = 2 * mu * (1 - L0)                  # leaf coefficient
    Y = Lm1[:, None] * W * U
    out += 12 * np.einsum('aj,aj->j', W2 * kk[:, None], c @ Y - L0[:, None] * W * (c2 @ (Lm1[:, None] * W)))
    t3 = (c2 * c) @ (X ** 3)
    out += 4 * np.einsum('aj,aj->j', Lm2[:, None] * W, U ** 3 - 3 * U * s + 2 * t3)
    out += 12 * np.einsum('bj,bj->j', Y, c @ Y)
    return out


def relu_map_edges4(m, C, K, k4):
    v = np.diag(C).copy(); k3 = np.diag(K).copy()
    L = ladder(m, v, k3, k4)
    mu, var, k3a, k4a = node_moments(L)
    c = C.copy(); np.fill_diagonal(c, 0.0)
    Kab = K.copy(); np.fill_diagonal(Kab, 0.0)
    T = edge_table(m, v, c, pmax=3)
    r, cc = (lambda x: x[:, None]), (lambda x: x[None, :])
    args = (r(k3), cc(k3), r(k4), cc(k4), Kab, Kab.T)
    E11 = pair_moment(T, 1, 1, *args); E21 = pair_moment(T, 2, 1, *args)
    E2 = 2 * L[2]; E3 = 6 * L[3]
    Ca = E11 - np.outer(mu, mu)
    Ka = E21 - np.outer(E2, mu) - 2 * mu[:, None] * Ca
    np.fill_diagonal(Ca, var); np.fill_diagonal(Ka, 0.0)
    Q, R = pair_cumulants4(T, mu, E2, E3, args)
    return mu, Ca, Ka, k3a, k4a, c, L, Q, R


def estimate_edge4(Ws, old=1, k4mode='full'):
    Ls, n, _ = Ws.shape
    W = Ws[0].astype(np.float64)
    m = np.zeros(n); C = W.T @ W; K = np.zeros((n, n)); k4 = np.zeros(n)
    out = []; prev = None
    for l in range(Ls):
        mu, Ca, Ka, k3a, k4a, c, L, Q, R = relu_map_edges4(m, C, K, k4)
        L0, Lm1 = L[0], L[-1]
        out.append(mu)
        if l + 1 == Ls:
            break
        Wn = Ws[l + 1].astype(np.float64)
        M = c * c * (L0 ** 2)[:, None] * Lm1[None, :]
        spec = (k3a, Ka - M, c, L0, Lm1)
        Kn = contract(*spec, Wn, Wn)
        if old and prev is not None:
            pspec, Wl = prev
            P = Wl @ (L0[:, None] * Wn)
            Kz = K.copy(); k3z = np.diag(K).copy(); np.fill_diagonal(Kz, 0.0)
            Kn = Kn + contract(*pspec, P, P) - contract(k3z * L0 ** 3, Kz * (L0 ** 2)[:, None] * L0[None, :],
                                                      None, L0, Lm1, Wn, Wn)
        prev = (spec, Wn)
        k4n = k4_next(Wn, k4a, Q, R, c, L, mu) if k4mode == 'full' else (Wn ** 4).T @ k4a
        m = mu @ Wn; C = Wn.T @ Ca @ Wn; K = Kn; k4 = k4n
    return np.array(out)


# ------------------------------------------------- v3: carry the order-parameter (Q) edge belief ----
def relu_map_v3(m, C, K, k4, Qz):
    v = np.diag(C).copy(); k3 = np.diag(K).copy()
    L = ladder(m, v, k3, k4)
    mu, var, k3a, k4a = node_moments(L)
    c = C.copy(); np.fill_diagonal(c, 0.0)
    Kab = K.copy(); np.fill_diagonal(Kab, 0.0)
    Qo = Qz.copy(); np.fill_diagonal(Qo, 0.0)
    T = edge_table(m, v, c, pmax=3)
    r, cc = (lambda x: x[:, None]), (lambda x: x[None, :])
    Z = np.zeros_like(c)
    args = (r(k3), cc(k3), r(k4), cc(k4), Kab, Kab.T, Qo, Z, Z)
    E11 = pair_moment(T, 1, 1, *args); E21 = pair_moment(T, 2, 1, *args)
    E2 = 2 * L[2]; E3 = 6 * L[3]
    Ca = E11 - np.outer(mu, mu)
    Ka = E21 - np.outer(E2, mu) - 2 * mu[:, None] * Ca
    np.fill_diagonal(Ca, var); np.fill_diagonal(Ka, 0.0)
    Q, R = pair_cumulants4(T, mu, E2, E3, args)
    return mu, Ca, Ka, k3a, k4a, c, L, Q, R


def estimate_v3(Ws, old=1, qgen=True):
    Ls, n, _ = Ws.shape
    W = Ws[0].astype(np.float64)
    m = np.zeros(n); C = W.T @ W; K = np.zeros((n, n)); k4 = np.zeros(n); Qz = np.zeros((n, n))
    out = []; prev = None
    for l in range(Ls):
        mu, Ca, Ka, k3a, k4a, c, L, Q, R = relu_map_v3(m, C, K, k4, Qz)
        L0, Lm1 = L[0], L[-1]
        out.append(mu)
        if l + 1 == Ls:
            break
        Wn = Ws[l + 1].astype(np.float64)
        M = c * c * (L0 ** 2)[:, None] * Lm1[None, :]
        spec = (k3a, Ka - M, c, L0, Lm1)
        Kn = contract(*spec, Wn, Wn)
        if old and prev is not None:
            pspec, Wl = prev
            P = Wl @ (L0[:, None] * Wn)
            Kz = K.copy(); k3z = np.diag(K).copy(); np.fill_diagonal(Kz, 0.0)
            Kn = Kn + contract(*pspec, P, P) - contract(k3z * L0 ** 3, Kz * (L0 ** 2)[:, None] * L0[None, :],
                                                      None, L0, Lm1, Wn, Wn)
        prev = (spec, Wn)
        W2 = Wn * Wn
        Qn = W2.T @ ((Q + np.diag(k4a)) @ W2)
        if qgen:
            X = L0[:, None] * Wn; U = c @ X; s = (c * c) @ (X * X)
            h4 = 2 * (L0 * (1 - L0) - mu * Lm1)
            H = (W2 * h4[:, None]).T @ (U * U - s)
            Qn = Qn + H + H.T
        k4n = k4_next(Wn, k4a, Q, R, c, L, mu)
        np.fill_diagonal(Qn, 0.0)
        m = mu @ Wn; C = Wn.T @ Ca @ Wn; K = Kn; k4 = k4n; Qz = Qn
    return np.array(out)


# ------------------------------------- v4: order-parameter fluctuation as a scale mixture (smooth) --
_QX = np.array([-np.sqrt(3.0), 0.0, np.sqrt(3.0)]); _QW = np.array([1 / 6, 2 / 3, 1 / 6])


def node_mixture_ladder(m, v, k3, k4, qmax=0.3):
    """Node belief with kappa4 carried as a variance scale mixture (q = kappa4 / 3v^2, 3-point GH),
    the residual kappa4 (outside [0, qmax]) by Edgeworth."""
    q = np.clip(k4 / (3 * v * v), 0.0, qmax)
    k4res = k4 - 3 * q * v * v
    Lm = None
    for x, w in zip(_QX, _QW):
        s = 1 + x * np.sqrt(q)
        Lx = ladder(m, v * s, k3 * s ** 1.5, k4res * s * s)
        Lm = {k: w * val for k, val in Lx.items()} if Lm is None else {k: Lm[k] + w * Lx[k] for k in Lm}
    return Lm


def relu_map_v4(m, C, K, k4, Qz, need4=True):
    """Pairs: E over a per-pair shared variance modulation (1 + eta), Var eta = q_ab = Q_ab/(v_a v_b),
    3-point Gauss-Hermite; inside each node the exact bivariate-Gaussian edge + node k3 + K Edgeworth."""
    v = np.diag(C).copy(); k3 = np.diag(K).copy()
    L = node_mixture_ladder(m, v, k3, k4)
    mu, var, k3a, k4a = node_moments(L)
    c = C.copy(); np.fill_diagonal(c, 0.0)
    Kab = K.copy(); np.fill_diagonal(Kab, 0.0)
    q = Qz / np.outer(v, v); np.fill_diagonal(q, 0.0)
    q = np.clip(q, 0.0, 0.3)
    r, cc = (lambda x: x[:, None]), (lambda x: x[None, :])
    keys = [(1, 1), (2, 1)] + ([(2, 2), (3, 1), (1, 2)] if need4 else [])
    Em = {k: 0.0 for k in keys}
    z0 = np.zeros_like(c)
    for x, w in zip(_QX, _QW):
        s = 1 + x * np.sqrt(q)
        va = v[:, None] * s; vb = v[None, :] * s
        T = edge_table_pairs(m, va, vb, c * s, pmax=3 if need4 else 2)
        args = (r(k3) * s ** 1.5, cc(k3) * s ** 1.5, z0, z0, Kab * s ** 1.5, Kab.T * s ** 1.5)
        for k in keys:
            Em[k] = Em[k] + w * pair_moment(T, k[0], k[1], *args)
        # mixture marginals of the same per-pair model (node k3 scaled with the variance)
        for side, mm, vv, kk in (('a', m[:, None], va, r(k3) * s ** 1.5), ('b', m[None, :], vb, cc(k3) * s ** 1.5)):
            G = uni_G(mm + 0 * vv, vv, -6, 3)
            for p in (1, 2, 3):
                val = factorial(p) * (G[p] + kk / 6 * G[p - 3])
                key = (p, 0) if side == 'a' else (0, p)
                Em[key] = Em.get(key, 0.0) + w * val
    ma, mb = Em[(1, 0)], Em[(0, 1)]
    Ca = Em[(1, 1)] - ma * mb
    if need4:
        Ea2 = Em[(2, 0)]
        Ka = Em[(2, 1)] - Ea2 * mb - 2 * ma * Ca
        from math import comb
        raw = lambda i, j: (1.0 if (i, j) == (0, 0) else Em[(i, j)] if (i, j) in Em else None)
        def raw2(i, j):
            if j == 0 and i == 0: return 1.0
            if (i, j) in Em: return Em[(i, j)]
            raise KeyError((i, j))
        def central(p, qq):
            return sum(comb(p, i) * comb(qq, j) * (-ma) ** (p - i) * (-mb) ** (qq - j) * raw2(i, j)
                       for i in range(p + 1) for j in range(qq + 1))
        c11 = central(1, 1); c20 = central(2, 0); c02 = central(0, 2)
        R = central(3, 1) - 3 * c20 * c11
        Q = central(2, 2) - c20 * c02 - 2 * c11 ** 2
        np.fill_diagonal(R, 0.0); np.fill_diagonal(Q, 0.0)
    else:
        Ka = Em[(2, 1)] - np.outer(2 * L[2], mu) - 2 * mu[:, None] * Ca
        Q = R = None
    np.fill_diagonal(Ca, var); np.fill_diagonal(Ka, 0.0)
    return mu, Ca, Ka, k3a, k4a, c, L, Q, R


def edge_table_pairs(m, va, vb, c, pmax=2, kmin=-4):
    """edge_table with per-pair variances va (n,n), vb (n,n)."""
    from math import comb
    ma, mb = m[:, None] + 0 * va, m[None, :] + 0 * vb
    rho = np.clip(c / np.sqrt(va * vb), -0.999999, 0.999999)
    c = rho * np.sqrt(va * vb)
    Ga = uni_G(ma, va, 2 * kmin - 1, pmax)
    Gb = uni_G(mb, vb, 2 * kmin - 1, pmax)
    vba = np.maximum(vb - c * c / va, 1e-12 * vb); mba = mb - c * ma / va
    vab = np.maximum(va - c * c / vb, 1e-12 * va); mab = ma - c * mb / vb
    gba = uni_G(mba, vba, 2 * kmin, pmax)
    gab = uni_G(mab, vab, 2 * kmin, pmax)
    T = {}
    for k in range(1, -kmin + 1):
        for j in range(kmin, pmax + 1):
            T[(-k, j)] = sum(comb(k - 1, i) * Ga[-k + i] * (-c / va) ** i * gba[j - i] for i in range(k))
            if j >= 0:
                T[(j, -k)] = sum(comb(k - 1, i) * Gb[-k + i] * (-c / vb) ** i * gab[j - i] for i in range(k))
    T[(0, 0)] = bvn_orthant(ma / np.sqrt(va), mb / np.sqrt(vb), rho)
    for qq in range(1, pmax + 1):
        T[(0, qq)] = (mb * T[(0, qq - 1)] + vb * T[(0, qq - 2)] + c * T[(-1, qq - 1)]) / qq
    for p in range(1, pmax + 1):
        for qq in range(0, pmax + 1):
            T[(p, qq)] = (ma * T[(p - 1, qq)] + va * T[(p - 2, qq)] + c * T[(p - 1, qq - 1)]) / p
    return T


def estimate_v4(Ws, old=1, qgen=True, verbose=False):
    Ls, n, _ = Ws.shape
    W = Ws[0].astype(np.float64)
    m = np.zeros(n); C = W.T @ W; K = np.zeros((n, n)); k4 = np.zeros(n); Qz = np.zeros((n, n))
    out = []; prev = None; info = []
    for l in range(Ls):
        mu, Ca, Ka, k3a, k4a, c, L, Q, R = relu_map_v4(m, C, K, k4, Qz, need4=(l + 1 < Ls))
        L0, Lm1 = L[0], L[-1]
        out.append(mu); info.append(dict(m=m, v=np.diag(C).copy(), k3=np.diag(K).copy(), k4=k4))
        if l + 1 == Ls:
            break
        Wn = Ws[l + 1].astype(np.float64)
        M = c * c * (L0 ** 2)[:, None] * Lm1[None, :]
        spec = (k3a, Ka - M, c, L0, Lm1)
        Kn = contract(*spec, Wn, Wn)
        if old and prev is not None:
            pspec, Wl = prev
            P = Wl @ (L0[:, None] * Wn)
            Kz = K.copy(); k3z = np.diag(K).copy(); np.fill_diagonal(Kz, 0.0)
            Kn = Kn + contract(*pspec, P, P) - contract(k3z * L0 ** 3, Kz * (L0 ** 2)[:, None] * L0[None, :],
                                                      None, L0, Lm1, Wn, Wn)
        prev = (spec, Wn)
        W2 = Wn * Wn
        Qn = W2.T @ ((Q + np.diag(k4a)) @ W2)
        if qgen:
            X = L0[:, None] * Wn; U = c @ X; s = (c * c) @ (X * X)
            h4 = 2 * (L0 * (1 - L0) - mu * Lm1)
            H = (W2 * h4[:, None]).T @ (U * U - s)
            Qn = Qn + H + H.T
        k4n = k4_next(Wn, k4a, Q, R, c, L, mu)
        np.fill_diagonal(Qn, 0.0)
        m = mu @ Wn; C = Wn.T @ Ca @ Wn; K = Kn; k4 = k4n; Qz = Qn
    return (np.array(out), info) if verbose else np.array(out)


def estimate_edge_info(Ws):
    """estimate_edge with per-layer node info (for diagnostics)."""
    Ls, n, _ = Ws.shape
    W = Ws[0].astype(np.float64)
    m = np.zeros(n); C = W.T @ W; K = np.zeros((n, n)); k4 = np.zeros(n)
    out = []; prev = None; info = []
    for l in range(Ls):
        mu, Ca, Ka, k3a, k4a, c, L0, Lm1 = relu_map_edges(m, C, K, k4)
        out.append(mu); info.append(dict(m=m, v=np.diag(C).copy(), k3=np.diag(K).copy(), k4=k4))
        if l + 1 == Ls:
            break
        Wn = Ws[l + 1].astype(np.float64)
        M = c * c * (L0 ** 2)[:, None] * Lm1[None, :]
        spec = (k3a, Ka - M, c, L0, Lm1)
        Kn = contract(*spec, Wn, Wn)
        if prev is not None:
            pspec, Wl = prev
            P = Wl @ (L0[:, None] * Wn)
            Kz = K.copy(); k3z = np.diag(K).copy(); np.fill_diagonal(Kz, 0.0)
            Kn = Kn + contract(*pspec, P, P) - contract(k3z * L0 ** 3, Kz * (L0 ** 2)[:, None] * L0[None, :],
                                                      None, L0, Lm1, Wn, Wn)
        prev = (spec, Wn)
        m = mu @ Wn; C = Wn.T @ Ca @ Wn; K = Kn; k4 = (Wn ** 4).T @ k4a
    return np.array(out), info
