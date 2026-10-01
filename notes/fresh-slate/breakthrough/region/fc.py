"""FC: first-order Wiener-chaos chain (prototype, numpy) — the sufficient construction tested at n = 1024.

Principle (REPORT.md Part 2): with fresh Gaussian weights each layer, the non-Gaussian structure that the final means
need is, to first order, the second Wiener chaos created at every ReLU and carried linearly (gated by Phi) to later
layers. A birth at layer s is the second-chaos component 1/2 w2_r :g_r^2: of a_{s,r} (w2 = E relu'' = phi/sigma); its
joint third cumulant with the first-chaos (Gaussian) part of z_l is fixed by two n x n "two-time" objects
    Y_s(l) = Cov(g_s, z_l) = S_s diag(P_s) Z_s(l)          (cross-covariance, Stein's lemma; exact for Gaussian g)
    Z_s(l) = W_{s+1} diag(P_{s+1}) W_{s+2} ... diag(P_{l-1}) W_l     (gated propagator of the chaos)
and D21(z_l)_ab = kappa3(z_a, z_a, z_b) = sum_s sum_r w2_{s,r} [ Y_ra^2 Z_rb + 2 Y_ra Z_ra Y_rb ].
Optional slice correction ('slices'): the exact Gaussian (2,1)/(3) slices of each birth replace the second-chaos ones,
carried as a slice-supported source with legs (Z, Z, Delta Z).

The covariance arrow uses the exact bivariate-Gaussian ReLU covariance plus the first-order bivariate Edgeworth term
G21 = sym(E[delta(z_a) 1(z_b > 0)] D21_ab); the per-neuron readout adds D3/6 E relu''' (D3 = diag D21) and the variance
gets the matching first-order term.

window w: sources older than w layers are dropped (w = None: all ages);  rank k: old sources (age > young) are carried
in a fixed row basis Q_s (n x k) = top-k left singular vectors of Z_s at the age they turn old.
"""
import numpy as np
from scipy.special import ndtr

import gclose as g

SQ2PI = np.sqrt(2 * np.pi)


def slice_exact(mu, S, P, m, nq=40):
    """Exact Gaussian slices of kappa3(a): K21_ik = kappa3(a_i, a_i, a_k) (i != k) and k3_i (diagonal), via 1-D
    Gauss-Hermite in z_i with the conditional ReLU mean of z_k | z_i in closed form."""
    v = np.diag(S).copy(); s = np.sqrt(v)
    x, wq = np.polynomial.hermite_e.hermegauss(nq); wq = wq / wq.sum()
    n = len(mu)
    K21 = np.zeros((n, n))
    k3 = np.zeros(n)
    cs = S.copy(); np.fill_diagonal(cs, 0)
    beta = cs / v[:, None]                       # regression coefficient of z_k on z_i (row i)
    condv = np.maximum(v[None, :] - cs * beta, 1e-12 * v[None, :])
    cs_ = np.sqrt(condv)
    for xq, wt in zip(x, wq):
        zi = mu + s * xq                         # (n,)
        fi = np.maximum(zi, 0) - m               # centred relu of z_i at this node
        cm = mu[None, :] + beta * (zi - mu)[:, None]          # conditional mean of z_k, rows i
        al = cm / cs_
        Ek = cm * ndtr(al) + cs_ * np.exp(-0.5 * al * al) / SQ2PI - m[None, :]   # E[f(z_k) - m_k | z_i]
        K21 += wt * (fi * fi)[:, None] * Ek
        k3 += wt * fi ** 3
    np.fill_diagonal(K21, 0)
    return K21, k3


def run(W, window=None, slices=False, young=None, rank=None, edge=True, readout=True, dtype=np.float32, trace=None, k4f=None, k4use=True):
    L, n, _ = W.shape
    W64 = W.astype(np.float64)
    Wf = W.astype(dtype)
    sources = []          # dicts: s, w2, Y, Z, (Delta), (Q)
    outs = []
    D21 = None
    mu = np.zeros(n); S = W64[0].T @ W64[0]
    for l in range(L):
        if l > 0:
            mu = m @ W64[l]; S = W64[l].T @ C @ W64[l]
        v = np.maximum(np.diag(S), 1e-300); s = np.sqrt(v); al = mu / s
        P = ndtr(al); p = np.exp(-0.5 * al * al) / SQ2PI
        mG = mu * P + s * p
        secG = (mu * mu + v) * P + mu * s * p
        C = g.cov_exact(mu, S, P)
        m = mG.copy(); var = secG - mG * mG
        if D21 is not None:
            D3 = np.diag(D21).copy()
            w3 = -al * p / s ** 2
            if readout:
                dm = D3 / 6 * w3
                m = mG + dm
                var = var + D3 * p / (3 * s) - 2 * mG * dm
            if k4f is not None and l in k4f and k4use:
                K22z, K31z, k4z = k4f[l]
                w4 = (al * al - 1) * p / s ** 3
                if readout:
                    dm4 = k4z / 24 * w4
                    var = var + k4z / 24 * 2 * w3 - 2 * m * dm4
                    m = m + dm4
                if edge:
                    X4 = 0.25 * np.outer(w2_, w2_) * K22z if False else 0.25 * np.outer(p / s, p / s) * K22z
                    X4 = X4 + (1 / 6) * (K31z * (w3[:, None] * P[None, :]))
                    X4 = X4 + X4.T - 0.25 * np.outer(p / s, p / s) * K22z
                    C += X4
            if edge:
                cs = S.copy(); np.fill_diagonal(cs, 0)
                va = v[:, None]; condv = np.maximum(v[None, :] - cs * cs / va, 1e-12 * v[None, :])
                E21w = (p / s)[:, None] * ndtr((mu[None, :] - cs * mu[:, None] / va) / np.sqrt(condv))
                X = E21w * D21
                C += 0.5 * (X + X.T)
        np.fill_diagonal(C, np.maximum(var, 1e-12))
        outs.append(m)
        if trace is not None:
            trace.append(dict(mu=mu, S=S, D21=D21, m=m, C=C))
        if l == L - 1:
            break
        # ---- chaos sources: birth at layer l, transport all to layer l + 1
        Wn = Wf[l + 1]
        G = (P.astype(dtype)[:, None] * Wn)                      # diag(P_l) W_{l+1}
        for src in sources:
            src["Z"] = src["Z"] @ G                               # (n, n) or, rank-k, R = Q^T Z (k, n)
        w2 = (p / s)
        new = dict(s=l, w2=w2.astype(dtype), Z=Wn.copy(), SP=(S * P[None, :]))
        if slices:
            K21x, k3x = slice_exact(mu, S, P, mG)
            # second-chaos (Wick) slices of the birth in a_l coordinates
            Cd = np.diag(S)
            Om = (w2[None, :] * (P[:, None] ** 2) * S.T ** 2) + 2 * (w2 * P * Cd)[:, None] * P[None, :] * S
            np.fill_diagonal(Om, 0)
            om3 = 3 * w2 * P ** 2 * Cd ** 2
            Dl = K21x - Om
            d3 = (k3x - om3)
            if D21 is not None and slices == 2:
                # first-order Edgeworth correction of the exact slices of kappa3(a_l) by the chain's own D21(z_l), minus
                # the part the gated transport of the old sources already carries (P_i^2 P_k D21_ik, P_i^3 D3_i)
                Dz = D21.copy(); D3z = np.diag(Dz).copy(); np.fill_diagonal(Dz, 0)
                w3 = -al * p / s ** 2
                M1 = mG; M2 = secG
                dM1 = D3z / 6 * w3; dM2 = D3z * w2 / 3; dM3 = D3z * P
                dk3 = dM3 - 3 * (dM1 * M2 + M1 * dM2) + 6 * M1 ** 2 * dM1
                Eik = C.copy(); np.fill_diagonal(Eik, 0); Eik = Eik + np.outer(M1, M1)       # E[a_i a_k], Gaussian
                Phi2 = np.outer(P, P) + S * np.outer(w2, w2) * 0 + (S - np.diag(np.diag(S))) * np.outer(p / s, p / s)
                dE2k = Dz * Phi2 + Dz.T * (M1[:, None] * w2[None, :]) + (D3z * w2 / 3)[:, None] * M1[None, :] \
                    + (M2[:, None] * w3[None, :]) * D3z[None, :] / 6
                dEik = 0.5 * Dz * (w2[:, None] * P[None, :]) + 0.5 * Dz.T * (P[:, None] * w2[None, :]) \
                    + (D3z * w3)[:, None] * M1[None, :] / 6 + M1[:, None] * (D3z * w3)[None, :] / 6
                dk21 = (dE2k - 2 * (dEik * M1[:, None] + Eik * dM1[:, None]) - (dM2[:, None] * M1[None, :] + M2[:, None] * dM1[None, :])
                        + 2 * (2 * (M1 * dM1)[:, None] * M1[None, :] + (M1 ** 2)[:, None] * dM1[None, :]))
                if k4f is not None and l in k4f:
                    K22z, K31z, k4z = k4f[l]
                    K22z = K22z.copy(); np.fill_diagonal(K22z, 0); K31z = K31z.copy(); np.fill_diagonal(K31z, 0)
                    w4 = (al * al - 1) * p / s ** 3
                    e1 = k4z * w4 / 24; e2 = k4z * w3 / 12; e3 = k4z * w2 / 4
                    dE2k4 = K31z * (w2[:, None] * P[None, :]) / 3 + 0.5 * K22z * (P[:, None] * w2[None, :]) \
                        + K31z.T * (M1[:, None] * w3[None, :]) / 3 + (k4z * w3 / 12)[:, None] * M1[None, :] \
                        + M2[:, None] * (k4z * w4 / 24)[None, :]
                    dEik4 = K31z * (w3[:, None] * P[None, :]) / 6 + K31z.T * (P[:, None] * w3[None, :]) / 6 \
                        + 0.25 * K22z * np.outer(w2, w2) + (k4z * w4 / 24)[:, None] * M1[None, :] + M1[:, None] * (k4z * w4 / 24)[None, :]
                    dm4 = e1; dM24 = e2; dM34 = e3
                    dk21 = dk21 + (dE2k4 - 2 * (dEik4 * M1[:, None] + Eik * dm4[:, None]) - (dM24[:, None] * M1[None, :] + M2[:, None] * dm4[None, :])
                                   + 2 * (2 * (M1 * dm4)[:, None] * M1[None, :] + (M1 ** 2)[:, None] * dm4[None, :]))
                    dk3 = dk3 + dM34 - 3 * (dm4 * M2 + M1 * dM24) + 6 * M1 ** 2 * dm4
                np.fill_diagonal(dk21, 0)
                Dl = Dl + dk21 - (P[:, None] ** 2) * P[None, :] * Dz
                d3 = d3 + dk3 - P ** 3 * D3z
                if trace is not None:
                    trace[-1]['sl21'] = K21x + dk21; trace[-1]['sl3'] = k3x + dk3; trace[-1]['sl21G'] = K21x.copy(); trace[-1]['sl3G'] = k3x
            np.fill_diagonal(Dl, d3 / 3)
            new["Delta"] = Dl.astype(dtype)
        sources.append(new)
        if window is not None:
            sources = [x for x in sources if l + 1 - x["s"] <= window]
        # ---- low-rank row basis for old sources: Z_s ~ Q_s R_s, Y_s = S_s P_s Z_s ~ A_s R_s, A_s = S_s P_s Q_s
        if rank is not None and young is not None:
            for src in sources:
                age = l + 1 - src["s"]
                if age == young + 1 and "Q" not in src:
                    U, sv, Vt = np.linalg.svd(src["Z"].astype(np.float64), full_matrices=False)
                    Q = U[:, :rank]
                    src["Q"] = Q.astype(dtype)
                    src["A"] = (src["SP"] @ Q).astype(dtype)
                    src["Z"] = (Q.T @ src["Z"].astype(np.float64)).astype(dtype)
                    if "Delta" in src:
                        src["DQ"] = (src["Delta"].astype(np.float64) @ Q).astype(dtype)
        # ---- D21 of z_{l+1}
        D = np.zeros((n, n), dtype=np.float64)
        for src in sources:
            w = src["w2"]
            if "Q" in src:
                Zf = src["Q"] @ src["Z"]; Y = src["A"] @ src["Z"]
                T = (src["DQ"] @ src["Z"]) if "DQ" in src else None
            else:
                Zf = src["Z"]; Y = src["SP"].astype(dtype) @ Zf
                T = (src["Delta"] @ Zf) if "Delta" in src else None
            D += (((Y * Y) * w[:, None]).T @ Zf).astype(np.float64)
            D += (2 * (((Y * Zf) * w[:, None]).T @ Y)).astype(np.float64)
            if T is not None:
                D += (((Zf * Zf).T @ T) + 2 * ((Zf * T).T @ Zf)).astype(np.float64)
        D21 = D
    return np.stack(outs)


def predict(W, **kw):
    return run(W, **kw)
