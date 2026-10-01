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

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../fresh-slate/breakthrough/region"))
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


def relu_k4(mu, v):
    """Exact fourth cumulant of relu(z), z ~ N(mu, v) (truncated-normal moments)."""
    s = np.sqrt(v); al = mu / s
    P = ndtr(al); p = np.exp(-0.5 * al * al) / SQ2PI
    I = [P, p, P - al * p, (al * al + 2) * p, -(al ** 3 + 3 * al) * p + 3 * P]
    from math import comb
    M = [None] + [sum(comb(k, j) * mu ** (k - j) * s ** j * I[j] for j in range(k + 1)) for k in range(1, 5)]
    m = M[1]
    c2 = M[2] - m * m
    c4 = M[4] - 4 * m * M[3] + 6 * m * m * M[2] - 3 * m ** 4
    return c4 - 3 * c2 * c2


def atom_norm2(w, Y, Z, T):
    """Exact ||A_r||_F^2 of each atom r: A_r = w_r[(y_r o y_r) z_r^T + 2 (y_r o z_r) y_r^T] + (z_r o z_r) t_r^T + 2 (z_r o t_r) z_r^T."""
    Y = Y.astype(np.float64); Z = Z.astype(np.float64); w = w.astype(np.float64)[:, None]
    A = [w * Y * Y, 2 * w * Y * Z]; B = [Z, Y]
    if T is not None:
        T = T.astype(np.float64); A += [Z * Z, 2 * Z * T]; B += [T, Z]
    out = np.zeros(Y.shape[0])
    for i in range(len(A)):
        for j in range(i, len(A)):
            t = (A[i] * A[j]).sum(1) * (B[i] * B[j]).sum(1)
            out += t if i == j else 2 * t
    return out


def poisson_q(sc, k):
    """Inclusion probabilities q = min(1, c sc) with sum q = k (water-filling)."""
    sc = np.maximum(sc, 1e-300); lo, hi = 0.0, 1e300
    c = k / sc.sum()
    for _ in range(100):
        q = np.minimum(1.0, c * sc); tot = q.sum()
        if abs(tot - k) < 1e-6 * k:
            break
        free = q < 1
        c *= (k - (~free).sum()) / max((c * sc[free]).sum(), 1e-300)
    return np.minimum(1.0, c * sc)


def run(W, window=None, slices=False, young=None, rank=None, edge=True, readout=True, dtype=np.float32, trace=None, k4f=None, k4use=True, share=None, share_young=2, k4own=False, k4mf=False, prune=None, prune_slices=None,
        sample=None, sample_young=2, seed=0, diag=None, diag_fracs=(0.5, 0.25), imp='norm', diag_store=None, oldproj=None, projsrc='prop', merge_age=None, merge_rep='oldest'):
    rng = np.random.default_rng(seed)
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
            if "Zs" in src:                                       # sampled old source: transport only its k rows
                src["Zs"] = src["Zs"] @ G; src["Ys"] = src["Ys"] @ G
                if src.get("Ts") is not None:
                    src["Ts"] = src["Ts"] @ G
                continue
            src["Z"] = src["Z"] @ G                               # (n, n) or, rank-k, R = Q^T Z (k, n)
        w2 = (p / s)
        new = dict(s=l, w2=w2.astype(dtype), Z=Wn.copy(), SP=(S * P[None, :]), k4a=relu_k4(mu, v))
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
        if prune is not None:
            nk = max(1, int(round(prune * n)))
            score = w2 * np.diag(S) * np.sqrt((Wn.astype(np.float64) ** 2).sum(1))
            new["Rw"] = np.sort(np.argsort(-score)[:nk])
        if prune_slices is not None and "Delta" in new:
            nk = max(1, int(round(prune_slices * n)))
            score = np.sqrt((new["Delta"].astype(np.float64) ** 2).sum(1)) * (Wn.astype(np.float64) ** 2).sum(1)
            new["Rs"] = np.sort(np.argsort(-score)[:nk])
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
        # ---- shared Oseledets subspace for old content: U_{l+1} = top-q right singular subspace of the oldest
        #      propagator; old sources (age > share_young) keep only the components of their legs inside it (dynamic)
        if share is not None:
            Zold = sources[0]["Z"].astype(np.float64)
            _, _, Vt = np.linalg.svd(Zold, full_matrices=False)
            U = Vt[:share].T.astype(dtype)                          # (n, q) basis of layer l+1 space
            for src in sources:
                if l + 1 - src["s"] > share_young:
                    src["Z"] = (src["Z"] @ U) @ U.T
                    if "Yp" not in src:
                        src["Yp"] = (src["SP"].astype(dtype) @ src["Z"])
                    else:
                        src["Yp"] = ((src["Yp"] @ G) @ U) @ U.T
        # ---- unbiased Poisson resampling of old atoms (Kyng-Sachdeva analogue): at age sample_young + 1 keep row r with
        #      probability q_r = min(1, c ||A_r||_F) (sum q = f n), reweight by 1/q_r, then transport only kept rows
        if sample is not None or diag is not None:
            for src in sources:
                age = l + 1 - src["s"]
                if age == sample_young + 1 and "Zs" not in src and "qd" not in src:
                    Zf = src["Z"]; Y = src["SP"].astype(dtype) @ Zf
                    T = (src["Delta"] @ Zf) if "Delta" in src else None
                    a2 = atom_norm2(src["w2"], Y, Zf, T)
                    sc = np.sqrt(a2) if imp == 'norm' else np.ones_like(a2)
                    if sample is not None:
                        q = poisson_q(sc, sample * n)
                        keep = rng.random(n) < q
                        c = (1.0 / q[keep]).astype(dtype)
                        src["Zs"] = Zf[keep]; src["Ys"] = Y[keep]; src["Ts"] = T[keep] if T is not None else None
                        src["c"] = c; src["ws"] = src["w2"][keep]; src["k"] = int(keep.sum())
                        del src["Z"]
                    else:
                        src["qd"] = {f: poisson_q(sc, f * n) for f in diag_fracs}
        # ---- D21 of z_{l+1}
        D21_prev = D21
        D = np.zeros((n, n), dtype=np.float64)
        Dold = np.zeros((n, n)) if diag is not None else None
        drec = dict(l=l + 1, E=0.0, sumA=0.0, var={f: 0.0 for f in diag_fracs}, varopt={f: [] for f in diag_fracs}, gsrc=[], nold=0)
        snaps = []
        for src in sources:
            if merge_age is not None:
                snaps.append(D.copy())
            w = src["w2"]
            if "Zs" in src:
                c = src["c"]; Y, Zr, Ts, wr = src["Ys"], src["Zs"], src["Ts"], src["ws"] * c
                D += (((Y * Y) * wr[:, None]).T @ Zr).astype(np.float64)
                D += (2 * (((Y * Zr) * wr[:, None]).T @ Y)).astype(np.float64)
                if Ts is not None:
                    D += ((((Zr * Zr) * c[:, None]).T @ Ts) + 2 * (((Zr * Ts) * c[:, None]).T @ Zr)).astype(np.float64)
                continue
            if diag is not None and "qd" in src:
                Zf = src["Z"]; Y = src["SP"].astype(dtype) @ Zf; T = (src["Delta"] @ Zf) if "Delta" in src else None
                Ds = (((Y * Y) * w[:, None]).T @ Zf).astype(np.float64) + 2 * (((Y * Zf) * w[:, None]).T @ Y).astype(np.float64)
                if T is not None:
                    Ds += (((Zf * Zf).T @ T) + 2 * ((Zf * T).T @ Zf)).astype(np.float64)
                a2 = atom_norm2(w, Y, Zf, T)
                if diag_store is not None and l + 1 in diag_store:
                    drec.setdefault("Ds", []).append(Ds.astype(np.float32)); drec.setdefault("Es", []).append(float(a2.sum()))
                Dold += Ds; drec["E"] += a2.sum(); drec["sumA"] += np.sqrt(a2).sum(); drec["nold"] += n
                drec["gsrc"].append(float((Ds ** 2).sum() / a2.sum()))
                for f in diag_fracs:
                    q = src["qd"][f]; drec["var"][f] += float(((1 / q - 1) * a2).sum())
                    drec["varopt"][f].append(a2)
                D += Ds
                continue
            if "Q" in src:
                Zf = src["Q"] @ src["Z"]; Y = src["A"] @ src["Z"]
                T = (src["DQ"] @ src["Z"]) if "DQ" in src else None
            else:
                Zf = src["Z"]; Y = src["Yp"] if "Yp" in src else src["SP"].astype(dtype) @ Zf
                T = (src["Delta"] @ Zf) if "Delta" in src else None
            if "Rw" in src:
                R = src["Rw"]; Yr, Zr, wr = Y[R], Zf[R], w[R]
            else:
                Yr, Zr, wr = Y, Zf, w
            D += (((Yr * Yr) * wr[:, None]).T @ Zr).astype(np.float64)
            D += (2 * (((Yr * Zr) * wr[:, None]).T @ Yr)).astype(np.float64)
            if T is not None:
                if "Rs" in src:
                    R = src["Rs"]; Zs, Ts = Zf[R], T[R]
                else:
                    Zs, Ts = Zf, T
                D += (((Zs * Zs).T @ Ts) + 2 * ((Zs * Ts).T @ Zs)).astype(np.float64)
        if merge_age is not None:
            # oracle age hierarchy: the block of sources with age >= merge_age is replaced by its least-squares fit on the
            # contributions of one or two representative sources (oldest; oldest + youngest of the block)
            snaps.append(D.copy())
            Dsl = [snaps[i + 1] - snaps[i] for i in range(len(sources))]
            blk = [i for i, x in enumerate(sources) if l + 1 - x["s"] >= merge_age]
            if len(blk) >= 2:
                Db = sum(Dsl[i] for i in blk)
                reps = [blk[0]] if merge_rep == 'oldest' else [blk[0], blk[-1]]
                Rm = np.stack([Dsl[i].ravel() for i in reps], 1)
                coef, *_ = np.linalg.lstsq(Rm, Db.ravel(), rcond=None)
                Dfit = (Rm @ coef).reshape(n, n)
                if diag is not None:
                    drec["merge_keep"] = float(1 - ((Db - Dfit) ** 2).sum() / (Db ** 2).sum()); drec["merge_coef"] = coef.tolist()
                D = D - Db + Dfit
            del snaps, Dsl
        if oldproj is not None and diag is not None and drec["nold"] > 0:
            # oracle: keep only the b-leg components of the old memory along q directions (Perron/Oseledets of the oldest
            # propagator, or the top right singular vectors of D_old itself)
            if projsrc == 'prop':
                _, _, Vt = np.linalg.svd(sources[0]["Z"].astype(np.float64), full_matrices=False)
            else:
                _, _, Vt = np.linalg.svd(Dold, full_matrices=False)
            V = Vt[:oldproj].T
            Dp = (Dold @ V) @ V.T
            drec["proj_keep"] = float((Dp ** 2).sum() / (Dold ** 2).sum())
            D = D - Dold + Dp
        D21 = D
        if diag is not None:
            drec["Dold2"] = float((Dold ** 2).sum()); drec["D2"] = float((D ** 2).sum())
            for f in diag_fracs:
                if drec["varopt"][f]:
                    a2 = np.concatenate(drec["varopt"][f]); q = poisson_q(np.sqrt(a2), f * len(a2))
                    drec["varopt"][f] = float(((1 / q - 1) * a2).sum())
                else:
                    drec["varopt"][f] = 0.0
            if diag_store is not None and l + 1 in diag_store:
                Z0 = sources[0]["Z"].astype(np.float64); v = np.ones(n)
                for _ in range(60):
                    v = Z0.T @ (Z0 @ v); v /= np.linalg.norm(v)
                drec["perron"] = v; drec["mu_next"] = None
                drec["Dold"] = Dold.astype(np.float32); drec["D"] = D.astype(np.float32)
            diag.append(drec)
        if k4mf:
            import k4mf as KM
            prev = k4f.get(l) if k4f is not None else None
            k4z = prev[2] if prev is not None else np.zeros(n)
            c22z = prev[0][0] if prev is not None else np.zeros(n)
            k4n, c22n = KM.step(mu, S, D21_prev, k4z, c22z, W64[l + 1])
            if k4f is None:
                k4f = {}
            k4f[l + 1] = (np.broadcast_to(c22n[None, :], (n, n)).copy(), np.zeros((n, n)), k4n)
        if k4own:
            c22 = np.zeros(n); k4d = np.zeros(n)
            for src in sources:
                Zf = (src["Q"] @ src["Z"]) if "Q" in src else src["Z"]
                Z2 = (Zf.astype(np.float64)) ** 2
                nu = Z2.sum(1)
                c22 += (src["k4a"] * nu) @ Z2 / n
                k4d += src["k4a"] @ (Z2 * Z2)
            if k4f is None:
                k4f = {}
            k4f[l + 1] = (np.broadcast_to(c22[None, :], (n, n)).copy(), np.zeros((n, n)), k4d)
    return np.stack(outs)


def predict(W, **kw):
    return run(W, **kw)
