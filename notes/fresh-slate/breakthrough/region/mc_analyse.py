"""Analyse an n = 1024 Monte Carlo atlas (mc1024.py): the local non-Gaussian error of the exact Gaussian closure at
every layer, its size in the norms that the linear-response transfer (lr_probe.py) prices, and how much of it the
joint third/fourth slices explain.

All squared norms are noise-corrected by the cross product of the two independent halves: ||X||^2 ~ <X_A, X_B>.

usage: python mc_analyse.py data/mc1024_mlp0_N262144.npz [results/lr_mlp0.json]
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
import gclose as g  # noqa: E402

OFF = None


def offd(X):
    Y = X.copy(); np.fill_diagonal(Y, 0); return Y


def cross(XA, XB):
    return float(np.sum(XA * XB))


def central(D, h, l):
    """Central moments of half h at layer l from raw moments about the reference."""
    d = D[f"{h}_s1y"][l]
    E2 = D[f"{h}_s2y"][l].astype(np.float64)
    E21 = D[f"{h}_s21"][l].astype(np.float64)          # E[y_a^2 y_b]
    E31 = D[f"{h}_s31"][l].astype(np.float64)          # E[y_a^3 y_b]
    E22 = D[f"{h}_s22"][l].astype(np.float64)          # E[y_a^2 y_b^2]
    E3 = D[f"{h}_s3"][l]; E4 = D[f"{h}_s4"][l]
    e2 = np.diag(E2).copy()
    C = E2 - np.outer(d, d)
    var = np.diag(C).copy()
    # kappa(X,X,Y) = E[X^2Y] - 2E[XY]E[X] - E[X^2]E[Y] + 2E[X]^2E[Y]
    K21 = E21 - 2 * E2 * d[:, None] - e2[:, None] * d[None, :] + 2 * (d * d)[:, None] * d[None, :]
    # central third/fourth along the diagonal
    m3 = E3 - 3 * d * e2 + 2 * d ** 3
    m4 = E4 - 4 * d * E3 + 6 * d * d * e2 - 3 * d ** 4
    k4 = m4 - 3 * var ** 2
    # central mu31 = E[(X-a)^3 (Y-b)]
    a = d[:, None]; b = d[None, :]
    mu31 = (E31 - 3 * a * E21 + 3 * a * a * E2 - a ** 3 * b - b * E3[:, None] + 3 * a * b * e2[:, None]
            - 3 * a * a * b * a + a ** 3 * b)
    # the two a^3 b terms: -a^3 E[Y] + a^3 b with E[Y] = b cancel; -3a^2 b E[X] + ... handled: E[X] = a
    K31 = mu31 - 3 * var[:, None] * C
    # mu22 = E[(X-a)^2 (Y-b)^2]
    mu22 = (E22 - 2 * b * E21 - 2 * a * E21.T + b * b * e2[:, None] + a * a * e2[None, :] + 4 * a * b * E2
            - 2 * a * b * b * a - 2 * a * a * b * b + a * a * b * b)
    # simplify the scalar tail exactly: terms in a^2 b^2: -2 -2 +1 ... recompute robustly below
    mu22 = E22 - 2 * b * E21 - 2 * a * E21.T + b * b * e2[:, None] + a * a * e2[None, :] + 4 * a * b * E2 - 3 * a * a * b * b
    K22 = mu22 - np.outer(var, var) - 2 * C * C
    return dict(d=d, C=C, var=var, K21=K21, k3=m3, k4=k4, K31=K31, K22=K22)


def analyse(path, lr_path=None):
    D = np.load(path)
    L = D["A_s1y"].shape[0]
    mu_ref = D["mu_ref"]; m_ref = D["m_ref"]
    lr = json.load(open(lr_path)) if lr_path and os.path.exists(lr_path) else None
    rows = []
    for l in range(L):
        per = {}
        for h in "AB":
            cm = central(D, h, l)
            mu = mu_ref[l] + cm["d"]
            S = cm["C"]
            mG, CG, P, p, s = g.closure_C(mu, S, "exact")
            t1 = D[f"{h}_s1t"][l]
            Ca = D[f"{h}_s2t"][l].astype(np.float64) - np.outer(t1, t1)
            mt = m_ref[l] + t1
            dC = Ca - CG
            al = mu / s
            w2 = p / s                                   # E f''  (gaussian)
            w3 = -al * p / s ** 2                        # E f''' = -alpha phi / sigma^2
            w4 = (al * al - 1) * p / s ** 3              # E f''''
            # first-order bivariate Edgeworth, factorised Gaussian weights (rho -> 0 weights)
            F21 = 0.5 * (w2[:, None] * P[None, :] * cm["K21"] + (w2[:, None] * P[None, :] * cm["K21"]).T)
            F31 = (1 / 6) * (w3[:, None] * P[None, :] * cm["K31"] + (w3[:, None] * P[None, :] * cm["K31"]).T)
            F22 = 0.25 * w2[:, None] * w2[None, :] * cm["K22"]
            # local mean error of the Gaussian readout and its first-order Edgeworth prediction
            dm = m_ref[l] - mG
            dm_pred = cm["k3"] / 6 * w3 + cm["k4"] / 24 * w4
            # exact bivariate-Gaussian weights for the (2,1) and (2,2) terms: E[delta(z_a) 1(z_b>0)], E[delta(z_a) delta(z_b)]
            v = np.diag(S).copy(); va = v[:, None]; vb = v[None, :]
            cs = S.copy(); np.fill_diagonal(cs, 0)
            condv = np.maximum(vb - cs * cs / va, 1e-12 * vb)
            E21w = (p / s)[:, None] * g.ndtr((mu[None, :] - cs * mu[:, None] / va) / np.sqrt(condv))
            det = np.maximum(va * vb - cs * cs, 1e-12 * va * vb)
            Q = (mu[:, None] ** 2 * vb + mu[None, :] ** 2 * va - 2 * cs * mu[:, None] * mu[None, :]) / det
            E22w = np.exp(-0.5 * Q) / (2 * np.pi * np.sqrt(det))
            G21 = 0.5 * (E21w * cm["K21"] + (E21w * cm["K21"]).T)
            G22 = 0.25 * E22w * cm["K22"]
            per[h] = dict(dC=dC, dCo=offd(dC), dvar=np.diag(dC).copy(), F21=offd(F21), F31=offd(F31), F22=offd(F22), G21=offd(G21), G22=offd(G22),
                          Coff=offd(Ca), dm=dm, dm_pred=dm_pred, mt=mt, mG=mG, K21=offd(cm["K21"]), var=cm["var"],
                          P=P)
        A, B = per["A"], per["B"]
        n = len(A["dm"])
        fro_dCo = cross(A["dCo"], B["dCo"])
        fro_Coff = cross(A["Coff"], B["Coff"])
        r = dict(layer=l, fro2_Coff=fro_Coff, fro2_dCoff=fro_dCo,
                 rms2_dvar=cross(A["dvar"], B["dvar"]) / n,
                 rms2_dm=cross(A["dm"], B["dm"]) / n,
                 fro2_K21=cross(A["K21"], B["K21"]),
                 fro2_F21=cross(A["F21"], B["F21"]), fro2_F31=cross(A["F31"], B["F31"]), fro2_F22=cross(A["F22"], B["F22"]))
        # coherent part of dC_off along C_off
        u = (A["Coff"] + B["Coff"]); u /= np.linalg.norm(u)
        r["coh_dCoff"] = 0.5 * (np.sum(A["dCo"] * u) + np.sum(B["dCo"] * u))
        # residual after first-order slices (fit-free, coefficient 1): cross-product of residuals
        for name, terms in [("21", ["F21"]), ("21+31", ["F21", "F31"]), ("21+31+22", ["F21", "F31", "F22"]), ("G21", ["G21"]), ("G21+31+G22", ["G21", "F31", "G22"])]:
            RA = A["dCo"] - sum(A[t] for t in terms); RB = B["dCo"] - sum(B[t] for t in terms)
            r["fro2_res_" + name] = cross(RA, RB)
        r["rms2_dm_res"] = cross(A["dm"] - A["dm_pred"], B["dm"] - B["dm_pred"]) / n
        if lr is not None and l < L - 1:
            r["mse_off"] = lr["off"][l] * fro_dCo
            r["mse_diag"] = lr["diag"][l] * r["rms2_dvar"]
            r["mse_mean"] = lr["mean"][l] * r["rms2_dm"]
        rows.append(r)
        print(json.dumps({k: (round(v, 6) if isinstance(v, float) and abs(v) > 1e-3 else (f"{v:.3e}" if isinstance(v, float) else v)) for k, v in r.items()}), flush=True)
    return rows


if __name__ == "__main__":
    rows = analyse(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
    out = os.path.join(HERE, "results", "mc_" + os.path.basename(sys.argv[1]).replace(".npz", ".json"))
    os.makedirs(os.path.dirname(out), exist_ok=True)
    json.dump(rows, open(out, "w"), indent=1)
