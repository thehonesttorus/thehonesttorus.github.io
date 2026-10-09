# Note XLIV section 5: tests T1-T3 on the production chain's cached state (noise-free via the two Monte Carlo halves).
#   python -I fw_test.py LOCDIR MC2DIR CD2DIR NET LADDER_PY
# LOCDIR: W_off{NET}.npy; MC2DIR: mc2_off{NET}_{full,h0,h1}.npz; CD2DIR: chaindump2_{NET}.npz;
# LADDER_PY: path to workbench/k3work/var_ladder.py (its helper definitions, KG, are reused).
# T1  R_s = Var_a(E_aa) / mean_(a!=b) E_ab^2 for E = C_chain(s) - C_true(s); the fresh-weight theorem predicts 2.
# T2  share of the omega-weighted D21 error explained by the dilation-fibre forms F1, F2, F3 (cross-fitted halves).
# T3  collective share c_k = 1 - <E0, Pk' E1 Pk'> / <E0, E1> of the off-diagonal error, k = 8, 32, 128; also for the
#     true non-Gaussian part g = C_true(s) - W_s KG(true state s-1) W_s^T.
import sys, numpy as np
_argv = list(sys.argv)
sys.argv = [_argv[5], ".", ".", "0", "1"]
exec(open(_argv[5]).read().split("cd = np.load")[0])    # ccoef, series, relu_var, KG, ndtr, ... (its header sets loc, net)
sys.argv = _argv
loc, mcd, cdd, net = _argv[1], _argv[2], _argv[3], int(_argv[4])   # parsed after the exec, which overwrites them
W = np.load(f"{loc}/W_off{net}.npy")
c2 = np.load(f"{cdd}/chaindump2_{net}.npz")
F = {h: np.load(f"{mcd}/mc2_off{net}_{h}.npz") for h in ("full", "h0", "h1")}
n = W.shape[1]
off = ~np.eye(n, dtype=bool)


def sym(A):
    return 0.5 * (A + A.T)


def chainC(s):
    C = sym(c2[f"C_off_{s}"].astype(np.float64)); np.fill_diagonal(C, c2[f"var_{s}"].astype(np.float64)); return C


def trueC(h, s):
    return sym(F[h]["cov"][s].astype(np.float64))


def ip(A, B):                                   # off-diagonal Frobenius inner product
    return float(np.sum(A[off] * B[off]))


def perp(E, U):                                 # P' E P' with P' = I - U U^T
    EU = E @ U; UE = U.T @ E
    return E - U @ UE - EU @ U.T + U @ (U.T @ EU) @ U.T


print(f"=== network {net}", flush=True)
print("T1: layer, R = Var_a(E_aa)/mean E_ab^2 (noise-free), rms relative variance error, common shift / rms", flush=True)
cache = {}
for s in range(2, 15):
    Cc = chainC(s)
    E = {h: Cc - trueC(h, s) for h in ("h0", "h1")}
    cache[s] = E
    d0, d1 = np.diag(E["h0"]), np.diag(E["h1"])
    m2 = np.mean(d0 * d1) - np.mean(d0) * np.mean(d1)
    o2 = ip(E["h0"], E["h1"]) / (n * (n - 1))
    v = np.mean(np.diag(trueC("full", s)))
    shift = 0.5 * (np.mean(d0) + np.mean(d1))
    print(f"  {s:2d}  R {m2 / o2:6.3f}   rms dv/v {np.sqrt(max(np.mean(d0 * d1), 0)) / v:.2e}   shift/rms "
          f"{shift / np.sqrt(max(m2, 1e-300)):+.3f}", flush=True)

print("T2: layer, explained share (F1+F2+F3, cross-fitted), single-form shares F1 / F2 / F3, noise-free corr(dD, chain D21)",
      flush=True)
for s in range(6, 15):
    mu = F["full"]["mu"][s].astype(np.float64); v = F["full"]["var"][s].astype(np.float64)
    C = trueC("full", s); sd = np.sqrt(v); al = mu / sd
    p0 = np.exp(-0.5 * al * al) / np.sqrt(2 * np.pi) / sd; Ph = ndtr(al)
    om = np.outer(p0, Ph); np.fill_diagonal(om, 0.0)
    forms = [np.outer(v + mu * mu, mu), C * mu[:, None], np.outer(mu * mu, mu)]
    X = [om * f for f in forms]
    Dc = c2[f"D21_{s}"].astype(np.float64)
    Y = {h: om * (F[h]["D21"][s].astype(np.float64) - Dc) for h in ("h0", "h1")}
    Enf = ip(Y["h0"], Y["h1"])
    G = np.array([[ip(a, b) for b in X] for a in X])
    b = {h: np.array([ip(a, Y[h]) for a in X]) for h in Y}
    share = float(b["h0"] @ np.linalg.solve(G, b["h1"])) / Enf
    single = [b["h0"][i] * b["h1"][i] / G[i, i] / Enf for i in range(3)]
    Xc = om * Dc
    corr = 0.5 * (ip(Xc, Y["h0"]) + ip(Xc, Y["h1"])) / np.sqrt(ip(Xc, Xc) * max(Enf, 1e-300))
    print(f"  {s:2d}  share {share:+.3f}   single {single[0]:+.3f} / {single[1]:+.3f} / {single[2]:+.3f}   "
          f"corr(error, chain D21) {corr:+.3f}", flush=True)

print("T3: layer, c_k of the chain's off-diagonal covariance error for k = 8 / 32 / 128 (random: 0.016 / 0.061 / 0.234);"
      " the same for the true non-Gaussian part g", flush=True)
for s in (12, 13, 14):
    Ct = trueC("full", s)
    lam, V = np.linalg.eigh(Ct); V = V[:, ::-1]
    E = {h: cache[s][h].copy() for h in cache[s]}
    for h in E:
        np.fill_diagonal(E[h], 0.0)
    g = {}
    for h in ("h0", "h1"):
        Kg = KG(F[h]["mu"][s - 1].astype(np.float64), trueC(h, s - 1))
        Ws = W[s].astype(np.float64)
        g[h] = trueC(h, s) - Ws @ Kg @ Ws.T; np.fill_diagonal(g[h], 0.0)
    e01, g01 = ip(E["h0"], E["h1"]), ip(g["h0"], g["h1"])
    ce, cg = [], []
    for k in (8, 32, 128):
        U = V[:, :k]
        ce.append(1.0 - ip(E["h0"], perp(E["h1"], U)) / e01)
        cg.append(1.0 - ip(g["h0"], perp(g["h1"], U)) / g01)
    print(f"  {s:2d}  error c_k {ce[0]:.3f} / {ce[1]:.3f} / {ce[2]:.3f}    true non-Gaussian part c_k "
          f"{cg[0]:.3f} / {cg[1]:.3f} / {cg[2]:.3f}    ||g||/||E|| {np.sqrt(max(g01, 0) / max(e01, 1e-300)):.2f}",
          flush=True)
