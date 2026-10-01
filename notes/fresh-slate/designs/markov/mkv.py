"""MKV(w): latent Gaussian Markov network with unary ReLU sites, depth junction tree of window w.

numpy float64 prototype (not metered).  See DESIGN.md §2.
  state: mu (n), C (n,n) Gaussian-field covariance, per-neuron k3, k4, window sources, old cumulants.
  w = number of layers a site's content is tracked exactly (w=1: only at the next layer); after that
  it is replaced by its Markov (per-neuron product) projection and transported diagonally.
"""
import numpy as np
from numpy.polynomial.hermite_e import hermeval
from scipy.special import ndtr

SQ2PI = np.sqrt(2 * np.pi)
_GL_X, _GL_W = np.polynomial.legendre.leggauss(80)


def pdf(t):
    return np.exp(-0.5 * t * t) / SQ2PI


def He(k, x):
    c = np.zeros(k + 1); c[k] = 1
    return hermeval(x, c)


# ------------------------------------------------------------------ per-neuron non-Gaussian ReLU moments
def relu_moments(mu, v, k3=None, k4=None):
    """E relu(z), E relu(z)^2, P(z>0) for z with mean mu, variance v, cumulants k3, k4 (Gram–Charlier,
    Edgeworth He6 term included)."""
    s = np.sqrt(v); t = mu / s; a = -t; ph = pdf(a)
    coef = {0: np.ones_like(mu)}
    if k3 is not None:
        g1 = k3 / s ** 3; g2 = k4 / v ** 2
        coef.update({3: g1 / 6, 4: g2 / 24, 6: g1 ** 2 / 72})

    def I0(k):
        return ndtr(t) if k == 0 else ph * He(k - 1, a)

    def I0m(k):  # with He_{-1} := 0 handled
        return 0.0 if k < 0 else I0(k)

    m1 = 0; m2 = 0; p = 0
    for k, c in coef.items():
        i0 = I0(k); i1 = I0(k + 1) + k * I0m(k - 1)
        i2 = I0(k + 2) + (2 * k + 1) * I0(k) + k * (k - 1) * I0m(k - 2)
        m1 = m1 + c * (mu * i0 + s * i1)
        m2 = m2 + c * (mu * mu * i0 + 2 * mu * s * i1 + v * i2)
        p = p + c * i0
    return m1, m2, p


# ------------------------------------------------------------------ site coefficients (Gaussian g)
def site_coeffs(mu, v, k3=None, k4=None):
    """Coefficients of the ReLU residual X(g) = relu(g) - E relu - E relu' (g - mu).
    g ~ N(mu, v) or, with k3/k4, the Gram-Charlier law with those cumulants (the site's actual
    marginal: the node potential of the Markov network).  E h^(k)(g) = s^-k int h f_k, where
    f_k = phi * sum_m c_m He_{m+k} (integration by parts on the Gram-Charlier density).
    Returns A3 = E X^3, A4 = k4(X), B21 = E (X^2)', B12 = E X'', B31 = E (X^3)', B22 = E (X^2)'',
    B13 = E X''' (derivatives in g)."""
    s = np.sqrt(v); t = mu / s
    cm = {0: np.ones_like(mu)}
    if k3 is not None:
        g1 = k3 / s ** 3; g2 = k4 / v ** 2
        cm.update({3: g1 / 6, 4: g2 / 24, 6: g1 ** 2 / 72})
    m1, _, p = relu_moments(mu, v, k3, k4) if k3 is not None else relu_moments(mu, v)
    u0 = -t
    lo, hi = -14.0, 14.0
    u0c = np.clip(u0, lo + 1e-9, hi - 1e-9)
    out = {k: 0.0 for k in ["X1", "X2", "A3", "X4", "B21", "B12", "B31", "B22", "B13"]}
    for (aa, bb) in [(np.full_like(u0c, lo), u0c), (u0c, np.full_like(u0c, hi))]:
        half = (bb - aa) / 2; mid = (bb + aa) / 2
        u = mid[:, None] + half[:, None] * _GL_X[None, :]
        wq = half[:, None] * _GL_W[None, :] * pdf(u)
        f = [sum(c[:, None] * He(m + k, u) for m, c in cm.items()) for k in range(4)]
        g = mu[:, None] + s[:, None] * u
        X = np.maximum(g, 0) - m1[:, None] - p[:, None] * s[:, None] * u
        X2 = X * X; X3 = X2 * X
        E = lambda h, k: (wq * h * f[k]).sum(1) / s[:] ** k
        out["X1"] += E(X, 0)
        out["X2"] += E(X2, 0); out["A3"] += E(X3, 0); out["X4"] += E(X2 * X2, 0)
        out["B21"] += E(X2, 1); out["B12"] += E(X, 2)
        out["B31"] += E(X3, 1); out["B22"] += E(X2, 2); out["B13"] += E(X, 3)
    out["A4"] = out["X4"] - 3 * out["X2"] ** 2
    return out


# ------------------------------------------------------------------ Gaussian bivariate ReLU covariance
def relu_cov_offdiag(mu, C, K=60):
    """Cov(relu z_c, relu z_d) for Gaussian (mu, C) by the Mehler series
    sum_k b_k^c b_k^d rho^k, b_k = s^k E[relu^(k)]/sqrt(k!)."""
    s = np.sqrt(np.diag(C)); t = mu / s
    R = C / np.outer(s, s); np.fill_diagonal(R, 0)
    out = np.zeros_like(C); Rk = np.ones_like(R); lf = 0.0
    for k in range(1, K + 1):
        Rk = Rk * R; lf += np.log(k)
        if k == 1:
            a = ndtr(t) * s
        else:
            a = He(k - 2, -t) * pdf(t) * s
        b = a * np.exp(-0.5 * lf)
        out += np.outer(b, b) * Rk
    return out


def single_site(P, K, c, Cm=None):
    P2 = P * P; P3 = P2 * P; K2 = K * K
    d3 = P3.T @ c["A3"] + 3 * (P2 * K).T @ c["B21"] + 3 * (P * K2).T @ c["B12"]
    d4 = (P2 * P2).T @ c["A4"] + 4 * (P3 * K).T @ c["B31"] + 6 * (P2 * K2).T @ c["B22"] \
        + 4 * (P * K2 * K).T @ c["B13"]
    if Cm is not None:
        # two-site trees (same source layer, one covariance edge c-d, c != d)
        Q1 = P2 * c["B21"][:, None]; Q2 = P * K * c["B12"][:, None]
        Co = Cm - np.diag(np.diag(Cm))
        CQ1, CQ2 = Co @ Q1, Co @ Q2
        d4 = d4 + 3 * (Q1 * CQ1).sum(0) + 12 * (Q2 * CQ1).sum(0) + 12 * (Q2 * CQ2).sum(0)
    return d3, d4


def mkv(W, w=4, ng_cov=True, K=60, return_all=False, gauss=False, var21=True, two_site=True, ng_site=False):
    """W: (L, n, n) float, x @ W convention. Returns per-layer means (L, n)."""
    W = np.asarray(W, dtype=np.float64)
    L, n, _ = W.shape
    mu = np.zeros(n); C = W[0].T @ W[0]; k3 = np.zeros(n); k4 = np.zeros(n)
    k3old = np.zeros(n); k4old = np.zeros(n)
    sources = []  # dicts: P, K, coef, age
    means = np.zeros((L, n)); diag = []
    for l in range(L):
        v = np.diag(C).copy()
        if ng_cov:
            m1, m2, p = relu_moments(mu, v, k3, k4)
        else:
            m1, m2, p = relu_moments(mu, v)
        means[l] = m1
        diag.append(dict(mu=mu.copy(), v=v, k3=k3.copy(), k4=k4.copy()))
        if l == L - 1:
            break
        Wn = W[l + 1]
        # field
        Ca = relu_cov_offdiag(mu, C, K); np.fill_diagonal(Ca, m2 - m1 ** 2)
        Cn = Wn.T @ Ca @ Wn
        if var21 and sources:
            # (2,1)-slice correction of Var z_{l+1}: sum_{c!=d} U_cj V_dj kappa(z_c,z_c,z_d), one site at a time
            s_ = np.sqrt(v); dens = pdf(mu / s_) / s_
            U = Wn * dens[:, None]; V = Wn * p[:, None]; UV = U * V
            dv = np.zeros(n); dC = np.zeros((n, n)); tdiag = np.zeros(n)
            for sd in sources:
                P_, K_, c = sd["P"], sd["K"], sd["coef"]
                PV, KV = P_ @ V, K_ @ V
                PPU, PKU, KKU = (P_ * P_) @ U, (P_ * K_) @ U, (K_ * K_) @ U
                # per-site kappa(z_c,z_c,z_c) for the c=d terms (they belong to Var a_c, already exact)
                td = c["A3"] @ (P_ ** 3) + 3 * c["B21"] @ (P_ * P_ * K_) + 3 * c["B12"] @ (P_ * K_ * K_)
                tdiag += td
                if var21 == "full":
                    M = (PPU * c["A3"][:, None]).T @ PV + (PPU * c["B21"][:, None]).T @ KV \
                        + 2 * (PKU * c["B21"][:, None]).T @ PV + 2 * (PKU * c["B12"][:, None]).T @ KV \
                        + (KKU * c["B12"][:, None]).T @ PV
                    dC += M
                else:
                    dv += c["A3"] @ (PPU * PV) + c["B21"] @ (PPU * KV + 2 * PKU * PV) \
                        + c["B12"] @ (2 * PKU * KV + KKU * PV)
            if var21 == "full":
                dC = dC + dC.T - (U.T @ (V * tdiag[:, None]) + (V * tdiag[:, None]).T @ U)
                Cn += 0.5 * dC
            else:
                dv -= ((U * V) * tdiag[:, None]).sum(0)
                Cn[np.diag_indices(n)] += dv
        mun = m1 @ Wn
        # sites at layer l and transport
        coef = site_coeffs(mu, v, k3, k4) if ng_site else site_coeffs(mu, v)
        newsrc = dict(P=Wn.copy(), K=(C * p[None, :]) @ Wn, coef=coef, age=1, Cm=C.copy() if two_site else None)
        for sdict in sources:
            sdict["P"] = (sdict["P"] * p[None, :]) @ Wn
            sdict["K"] = (sdict["K"] * p[None, :]) @ Wn
            sdict["age"] += 1
        sources.append(newsrc)
        # old content: diagonal (Markov) transport
        Wp = Wn * p[:, None]
        k3old = (Wp ** 3).T @ k3old; k4old = (Wp ** 4).T @ k4old
        if gauss:
            sources = []
        k3 = k3old.copy(); k4 = k4old.copy(); keep = []
        for sdict in sources:
            d3, d4 = single_site(sdict["P"], sdict["K"], sdict["coef"], sdict["Cm"])
            k3 += d3; k4 += d4
            if sdict["age"] >= w:
                k3old += d3; k4old += d4
            else:
                keep.append(sdict)
        sources = keep
        mu, C = mun, Cn
    return (means, diag) if return_all else means
