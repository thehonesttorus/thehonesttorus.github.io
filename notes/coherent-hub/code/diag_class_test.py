# Note XLV section 4: is the truth's non-package (3,1) content the transported single-site class sum_a W_ia^3 W_ja kappa4(h_a)?
#   python -I diag_class_test.py LOCDIR MC2DIR CD2DIR NET
# kappa_r(h_a), h = relu(z), from the true layer s-1 marginal of z_a (mean, variance, kappa3, kappa4; Edgeworth density on a
# 120-point Gauss-Hermite rule). Candidates for K31'[i,j] = kappa(z'_i, z'_i, z'_i, z'_j) at layer s:
#   X4  = sum_a W_ia^3 W_ja kappa4(h_a)                      (the (4) class of h: one site, four legs)
#   f31 = 12 v_i C'_ij                                       (dilation package)
# and the same for the (2,2) slice, X4_22 = sum_a W_ia^2 W_ja^2 kappa4(h_a), and for D21' (X3 = sum_a W_ia^2 W_ja kappa3(h_a)).
# Fits on neurons with true alpha > -2.5 (active by active, off-diagonal); explained energies noise-free across halves;
# chain error = truth - chain (chaindump2) regressed on each candidate.
import sys, numpy as np
loc, mcd, cdd, net = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
f64 = lambda a: np.asarray(a, dtype=np.float64)
sym = lambda A: 0.5 * (A + A.T)
W = np.load(f"{loc}/W_off{net}.npy").astype(np.float64)
c2 = np.load(f"{cdd}/chaindump2_{net}.npz")
T = {h: np.load(f"{mcd}/mc2_off{net}_{h}.npz") for h in ("full", "h0", "h1")}
n = W.shape[1]
off = ~np.eye(n, dtype=bool)
xg, wg = np.polynomial.hermite_e.hermegauss(120); wg = wg / wg.sum()
He3 = xg ** 3 - 3 * xg; He4 = xg ** 4 - 6 * xg ** 2 + 3; He6 = xg ** 6 - 15 * xg ** 4 + 45 * xg ** 2 - 15


def zd(A):
    A = A.copy(); np.fill_diagonal(A, 0.0); return A


def relu_cumulants(mu, v, k3, k4):
    S = np.sqrt(v); g3 = k3 / S ** 3; g4 = k4 / S ** 4
    dens = 1.0 + np.outer(g3, He3) / 6.0 + np.outer(g4, He4) / 24.0 + np.outer(g3 * g3, He6) / 72.0
    w = wg[None, :] * dens
    w = w / w.sum(axis=1, keepdims=True)
    h = np.maximum(mu[:, None] + S[:, None] * xg[None, :], 0.0)
    m = (w * h).sum(1); d = h - m[:, None]
    c2_, c3_, c4_ = (w * d ** 2).sum(1), (w * d ** 3).sum(1), (w * d ** 4).sum(1)
    return c2_, c3_, c4_ - 3 * c2_ ** 2


print(f"=== network {net}: single-site classes transported by the fresh rows, against the truth and the chain's error")
print("slice: R2 package alone -> package + X (noise-free), c_X;  X alone R2;  chain error ~ c X: c, share of error energy")
for s in (4, 6, 8, 10, 12, 13, 14):
    Ff = T["full"]
    mu0, v0, k30, k40 = (f64(Ff[k][s - 1]) for k in ("mu", "var", "k3", "k4"))
    _, kh3, kh4 = relu_cumulants(mu0, v0, k30, k40)
    Ws = W[s]; W2 = Ws * Ws
    X = {"K31": zd(((Ws * W2) * kh4[None, :]) @ Ws.T), "K22": zd((W2 * kh4[None, :]) @ W2.T),
         "D21": zd((W2 * kh3[None, :]) @ Ws.T)}
    mu1 = f64(Ff["mu"][s]); C1 = sym(f64(Ff["cov"][s])); v1 = np.diag(C1).copy(); C1o = zd(C1)
    act = (mu1 / np.sqrt(v1)) > -2.5
    mk = off & act[:, None] & act[None, :]
    forms = dict(K31=12.0 * v1[:, None] * C1o, K22=zd(4.0 * np.outer(v1, v1) + 8.0 * C1o * C1o),
                 D21=zd(2.0 * (2.0 * mu1[:, None] * C1 + mu1[None, :] * v1[:, None])))
    chain = dict(K31=zd(f64(c2[f"wk431_{s}"]).T), K22=zd(sym(f64(c2[f"wk4m_{s}"]))), D21=zd(f64(c2[f"D21_{s}"])))
    ip = lambda A, B: float(np.sum((A * B)[mk]))
    out = [f"  s={s:2d}"]
    for key in ("K31", "K22", "D21"):
        tr = {h: (zd(f64(T[h][key][s])) if key != "K22" else zd(sym(f64(T[h][key][s])))) for h in ("full", "h0", "h1")}
        f, Xk = forms[key], X[key]
        ET = ip(tr["h0"], tr["h1"])
        cp0 = ip(f, tr["full"]) / ip(f, f)
        r2p = 1 - ip(tr["h0"] - cp0 * f, tr["h1"] - cp0 * f) / ET
        A2 = np.array([[ip(f, f), ip(f, Xk)], [ip(f, Xk), ip(Xk, Xk)]]); b2 = np.array([ip(f, tr["full"]), ip(Xk, tr["full"])])
        cp, cx = np.linalg.solve(A2, b2)
        r2j = 1 - ip(tr["h0"] - cp * f - cx * Xk, tr["h1"] - cp * f - cx * Xk) / ET
        cx0 = ip(Xk, tr["full"]) / ip(Xk, Xk)
        r2x = 1 - ip(tr["h0"] - cx0 * Xk, tr["h1"] - cx0 * Xk) / ET
        e = {h: tr[h] - chain[key] for h in tr}
        EE = ip(e["h0"], e["h1"]); ce = ip(e["full"], Xk) / ip(Xk, Xk)
        she = 1 - ip(e["h0"] - ce * Xk, e["h1"] - ce * Xk) / EE
        # chain error with both package and X
        A3 = np.array([[ip(f, f), ip(f, Xk)], [ip(f, Xk), ip(Xk, Xk)]]); b3 = np.array([ip(f, e["full"]), ip(Xk, e["full"])])
        ep, ex = np.linalg.solve(A3, b3)
        she2 = 1 - ip(e["h0"] - ep * f - ex * Xk, e["h1"] - ep * f - ex * Xk) / EE
        out.append(f"    {key}: R2 pkg {100 * r2p:5.1f} -> +X {100 * r2j:5.1f} (c_p {cp:+.2f}, c_X {cx:+.2f}) | X alone {100 * r2x:5.1f} "
                   f"(c {cx0:+.2f}) | chain err: c {ce:+.2f} share {100 * she:5.1f}%, with pkg {100 * she2:5.1f}% (c_X {ex:+.2f})")
    print("\n".join(out), flush=True)
