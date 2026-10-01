"""v1: Mehler copula + one carried cumulant hyperedge (the pairwise (2,1) slice D_ac = kappa(z_a, z_a, z_c)).

The copula generates all joint structure from (T, R); the measured dominant defect of v0 is the (2,1) slice
(DESIGN.md section 4 / RESULTS.md). v1 carries D and corrects every bivariate expectation the step needs by the
first-order hyperedge (Edgeworth) term, evaluated in the same Mehler/hafnian basis:
  E[phi(z)] += (D - D_cop)_{ac}/2 E_cop[d_a^2 d_c phi] + (D - D_cop)_{ca}/2 E_cop[d_a d_c^2 phi].
"""
import numpy as np
from numpy.polynomial import polynomial as P
from copula import HE, FACT, K, PQ, upper_moments, positive_intervals, fleishman, latent_corr

FULL = upper_moments(-np.inf, 40)


def ext_profiles(m, s, a, dmax=K):
    """per neuron: h (F), F2, F3, F4 profiles; Hd (indicator 1{T>0}); dl (delta(T)); t (T); t2c ((T-m)^2)."""
    n = len(m)
    out = {k: np.zeros((n, dmax + 1)) for k in ['h', 'F2', 'F3', 'F4', 'Hd', 'dl', 't', 't2c']}
    for i in range(n):
        T = np.array([m[i] - s[i] * a[i, 1], s[i] * (a[i, 0] - 3 * a[i, 2]), s[i] * a[i, 1], s[i] * a[i, 2]])
        Tc = T.copy(); Tc[0] -= m[i]
        ivs = positive_intervals(T)
        kmax = 12 + dmax + 1
        Mint = np.zeros(kmax + 1)
        for lo, hi in ivs:
            Mint += upper_moments(lo, kmax) - upper_moments(hi, kmax)
        Tc_ = np.trim_zeros(T, 'b')
        r = np.roots(Tc_[::-1]) if len(Tc_) > 1 else np.array([])
        r = r[np.abs(r.imag) < 1e-9].real
        dT = P.polyder(T)
        Tr = [np.array([1.0])]
        for _ in range(4):
            Tr.append(P.polymul(Tr[-1], T))
        T2c = P.polymul(Tc, Tc)
        for d in range(dmax + 1):
            for k, key in zip(range(1, 5), ['h', 'F2', 'F3', 'F4']):
                c = P.polymul(Tr[k], HE[d]); out[key][i, d] = c @ Mint[:len(c)]
            out['Hd'][i, d] = HE[d] @ Mint[:len(HE[d])]
            out['dl'][i, d] = sum(P.polyval(x, HE[d]) * np.exp(-x * x / 2) / np.sqrt(2 * np.pi) / abs(P.polyval(x, dT)) for x in r)
            c = P.polymul(T, HE[d]); out['t'][i, d] = c @ FULL[:len(c)]
            c = P.polymul(T2c, HE[d]); out['t2c'][i, d] = c @ FULL[:len(c)]
    return out


def M(fa, fb, R, d0=1):
    out = np.zeros_like(R)
    Rp = np.ones_like(R) if d0 == 0 else R ** d0
    for d in range(d0, min(fa.shape[1], fb.shape[1])):
        out += np.outer(fa[:, d], fb[:, d]) * Rp / FACT[d]
        Rp = Rp * R
    return out


def layer_step(m, s, a, R, D, W, cfg):
    n = len(m)
    pr = ext_profiles(m, s, a)
    h, F2, F3, F4, Hd, dl = pr['h'], pr['F2'], pr['F3'], pr['F4'], pr['Hd'], pr['dl']
    mu = h[:, 0]
    e0 = np.eye(1, K + 1, 0)
    A2 = F2 - 2 * mu[:, None] * h + (mu ** 2)[:, None] * e0
    A3 = F3 - 3 * mu[:, None] * F2 + 3 * (mu ** 2)[:, None] * h - (mu ** 3)[:, None] * e0
    var = A2[:, 0]; k3a = A3[:, 0]
    k4a = F4[:, 0] - 4 * mu * F3[:, 0] + 6 * mu ** 2 * F2[:, 0] - 3 * mu ** 4 - 3 * var ** 2
    R0 = R.copy(); np.fill_diagonal(R0, 0.0)
    off = ~np.eye(n, dtype=bool)
    # hyperedge defect of the copula in the (2,1) slice of z
    Dcop = M(pr['t2c'], pr['t'], R0)
    Delta = (D - Dcop) * off if cfg.get('hyper', 1) else np.zeros((n, n))
    DeltaT = Delta.T
    # covariance of a
    Ca = M(h, h, R0) + 0.5 * Delta * M(dl, Hd, R0, 0) + 0.5 * DeltaT * M(Hd, dl, R0, 0)
    Ca[np.diag_indices(n)] = var
    # (2,1) slice of a
    hm = h - mu[:, None] * Hd
    K21 = M(A2, h, R0) + Delta * M(Hd - mu[:, None] * dl, Hd, R0, 0) + DeltaT * M(hm, dl, R0, 0)
    K21 = K21 * off
    mz = mu @ W
    Cz = W.T @ Ca @ W
    vz = np.diag(Cz).copy()
    W2, W3, W4 = W ** 2, W ** 3, W ** 4
    # paths (copula, distinct indices)
    X = {p: (R0 ** p) @ (W * h[:, p:p + 1]) / FACT[p] for p in range(1, PQ + 1)}
    t111 = np.zeros(W.shape[1])
    for p in range(1, PQ + 1):
        for q in range(1, PQ + 1):
            Z = (R0 ** (p + q)) @ (W2 * (h[:, p] * h[:, q])[:, None]) / (FACT[p] * FACT[q])
            t111 += np.sum(W * h[:, p + q][:, None] * (X[p] * X[q] - Z), 0)
    K21W = K21 @ W
    k3z = k3a @ W3 + 3 * np.sum(W2 * K21W, 0) + 3 * t111
    # fourth cumulant (copula sectors as v0)
    K22 = M(A2, A2, R0) - 2 * Ca ** 2 * off
    K31 = M(A3, h, R0) - 3 * var[:, None] * (Ca * off)
    t211 = np.zeros(W.shape[1])
    for p in range(1, PQ + 1):
        for q in range(1, PQ + 1):
            G = A2[:, p + q] - 2 * h[:, p] * h[:, q]
            Zbc = (R0 ** (p + q)) @ (W2 * (h[:, p] * h[:, q])[:, None]) / (FACT[p] * FACT[q])
            t211 += np.sum(W2 * G[:, None] * (X[p] * X[q] - Zbc), 0)
            Y = (R0 ** p) @ (W2 * A2[:, p:p + 1]) / FACT[p]
            t211 += 2 * np.sum(W * h[:, p + q][:, None] * Y * X[q], 0)
    k4z = k4a @ W4 + 3 * np.sum(W2 * ((K22 * off) @ W2), 0) + 4 * np.sum(W3 * ((K31 * off) @ W), 0) + 6 * t211
    # new (2,1) slice D'_{jk} = sum W_aj W_bj W_ck kappa(a_a, a_b, a_c)
    Dn = W2.T @ (K21W + k3a[:, None] * W)                       # a = b (incl. a = b = c)
    Dn += 2 * (W * K21W).T @ W                                   # a = c != b and b = c != a
    # all distinct, copula paths p = q = 1 (+ higher p, q)
    for p in range(1, PQ + 1):
        for q in range(1, PQ + 1):
            Xc = X[p] * X[q]                                     # centre c, both j-legs at the ends
            Dn += ((h[:, p + q][:, None] * Xc).T @ W)
            Dn += 2 * (W * h[:, p + q][:, None] * X[p]).T @ X[q]  # centre a (a j-leg), ends b (j) and c (k)
    Dn = Dn * (~np.eye(Dn.shape[0], dtype=bool)) + np.diag(k3z)
    sz = np.sqrt(vz)
    g1 = np.clip(k3z / sz ** 3, -1.5, 1.5); g2 = np.clip(k4z / sz ** 4, -1.0, 4.0)
    an = fleishman(g1, g2)
    t = np.stack([sz * an[:, 0], 2 * sz * an[:, 1], 6 * sz * an[:, 2]], 1)
    Rn = latent_corr(Cz, t)
    return mu, (mz, sz, an, Rn, Dn), dict(k3=k3z, k4=k4z, vz=vz, Cz=Cz, D=Dn)


def estimate(W, cfg=None, return_diag=False):
    cfg = cfg or {}
    L, n, _ = W.shape
    C1 = W[0].T @ W[0]
    s = np.sqrt(np.diag(C1)); R = C1 / np.outer(s, s)
    m = np.zeros(n); a = np.tile([1.0, 0.0, 0.0], (n, 1)); D = np.zeros((n, n))
    means, diags = [], []
    for l in range(L):
        if l < L - 1:
            mu, (m, s, a, R, D), dg = layer_step(m, s, a, R, D, W[l + 1].astype(float), cfg)
            diags.append(dg)
        else:
            mu = ext_profiles(m, s, a, dmax=0)['h'][:, 0]
        means.append(mu)
    return (np.array(means), diags) if return_diag else np.array(means)
