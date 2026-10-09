# The chaos-2 kappa4 diagonal against Monte Carlo truth (note XLI).   python chaos2_diag.py NET TAG
# Truncated at the second Wiener chaos of the input, a pre-activation is z_i = mu_i + L_i.x + (x^T H_i x - tr H_i) / 2,
# a generalized chi-square, whose fourth cumulant is 12 |H_i L_i|^2 + 3 tr H_i^4. The sources carry H: a source born at
# layer b contributes sum_m P_im w2_m l_m l_m^T (l_m = the first-chaos vector of z_(b,m), w2 = E f''). So
#   v_i = H_i L_i = sum_s sum_m X_(s,im) l_m,   X_s = P_s o At_s o w2_s,   At_s = A_s + P_s d(w1_b var_b)   (the full arm),
# and |v_i|^2 is the path class of the kappa4 diagonal for every pair of births. Three ways:
#   exact  : v_i from the first-chaos maps L_b (L_0 = W_0, L_(b+1) = W_(b+1) d(w1_b) L_b), all pairs;
#   within : 12 sum_s diag(X_s C_b X_s^T), births of one source only (one n^3 product per source-layer);
#   schur  : 12 diag(Y (C + eps)^-1 Y^T), Y = sum_s X_s At_s^T, all pairs through one solve (the identity
#            Cov(z_b) = Cov(z_b, z_l) Cov(z_l)^-1 Cov(z_l, z_b) of the first chaos when L_l is invertible).
# Also the chaos-3 star with the full arm, 4 sum c3 P At^3, against the V55 star 4 sum c3 P A^3.
# The chain runs the adopted fold system with every source dense (V21_NO_CONFINE=1), so all legs are available.
import sys, os, importlib.util, time, numpy as np, flopscope as flops
from whestbench import MLP
net, tag = int(sys.argv[1]), sys.argv[2]
LAYERS = [int(x) for x in os.environ.get("C2_LAYERS", "3,5,7,9,11,13,14,15").split(",")]
CYC_LAYERS = tuple(int(x) for x in os.environ.get("C2_CYC", "7,13").split(",") if x)
CYC_ROWS = list(range(0, 1024, 128))
EPS = [float(x) for x in os.environ.get("C2_EPS", "0,1e-4,1e-3,1e-2,3e-2,1e-1").split(",")]
_prod = {"V29_WARM_JOIN": "1", "V29_WARM_FB": "1", "V17_R_RES": "4", "V26_STRASSEN": "6", "V26_STRASSEN_MIN": "16",
         "V32_JOIN_SMM": "1", "V32_ROT_SMM": "1", "V32_JOIN_POST": "3", "V32_JP_C": "0.1", "V21_R_OLD": "320",
         "V24_R_OLD2": "192", "V18_R_FB": "2", "V24_AGE_OLD2": "8", "V33_K4Q": "3", "V33_K4Q_RANK": "4",
         "V34_OPT": "abcdefgj", "V33_SAT": "-2.5", "V35_SAT_ROUND": "64", "V35_SATC": "1", "V35_SKIP_JOIN_L2": "1",
         "V35_SB_MN": "8", "V35_CPRE_MN": "8", "V35_JN_MN": "8", "V52_FB_FOLD": "1", "V21_NO_CONFINE": "1"}
for _k, _v in _prod.items():
    os.environ.setdefault(_k, _v)
os.environ["V29_DUMP_LAYERS"] = ",".join(str(x) for x in range(16))
os.environ["V30_DUMP_LEGS"] = "1"
spec = importlib.util.spec_from_file_location("estv29", "est_v29.py"); mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
Wcol = np.load(f"../official/W_off{net}.npy"); mt = np.load(f"../official/truth_off{net}.npz")["m"].astype(np.float64)
mlp = MLP(width=1024, depth=16, weights=[np.ascontiguousarray(W.T).astype(np.float32) for W in Wcol])
t0 = time.time()
with flops.BudgetContext(flop_budget=2**45, wall_time_limit_s=3000.0, quiet=True) as _bc:
    out = np.asarray(mod.Estimator().predict(mlp, 2**41), dtype=np.float64)
print(f"net {net} {tag}: chain raw {np.mean((out[-1] - mt[-1]) ** 2):.5e}  ({time.time() - t0:.0f}s)", flush=True)
Wd = Wcol.astype(np.float64)
L_, n = Wd.shape[0], Wd.shape[1]

# per-layer Wick-stage quantities (first dict per layer that has them)
st = {}
for d in mod.DUMPS:
    l = d["layer"]
    for k in ("mu", "var", "C_off", "g4row", "w1", "w2", "D21", "D3"):
        if d.get(k) is not None and (k, l) not in st:
            st[(k, l)] = np.asarray(d[k], np.float64)
legs = {d["layer"]: d for d in mod.LEGS}
print("legs at layers", sorted(legs), " slots", {l: (None if legs[l]["A"] is None else legs[l]["A"].shape[0]) for l in sorted(legs)}, flush=True)

def phi(x):
    return np.exp(-0.5 * x * x) / np.sqrt(2 * np.pi)

# first-chaos maps L_b (z_b = W_b y_(b-1), L_0 = W_0); w1_b = E f'(z_b) as the chain uses it
Lc = [Wd[0]]
for b in range(1, L_):
    Lc.append(Wd[b] @ (st[("w1", b - 1)][:, None] * Lc[b - 1]))
mc = {h: np.load(f"mc2_off{net}_{h}.npz") for h in ("full", "h0", "h1")}

def stats(name, c, r_full, r_h0, r_h1):
    c = c - c.mean()
    S = np.mean((r_h0 - r_h0.mean()) * (r_h1 - r_h1.mean()))      # signal variance of the residual (split halves)
    cv = np.mean((r_full - r_full.mean()) * c)
    vc = np.mean(c * c)
    rho = cv / np.sqrt(max(np.var(r_full), 1e-300) * max(vc, 1e-300))
    rho_s = cv / np.sqrt(max(S, 1e-300) * max(vc, 1e-300)) if S > 0 else float("nan")
    b = cv / max(vc, 1e-300)
    frac = (cv * cv / max(vc, 1e-300)) / S if S > 0 else float("nan")   # share of the residual's signal variance removed
    return f"    {name:<16s} rms {np.sqrt(vc):.3e}  corr {rho:+.3f}  corr/signal {rho_s:+.3f}  slope {b:+.3f}  signal removed {frac:+.3f}"

for l in LAYERS:
    if l not in legs or legs[l]["A"] is None:
        print(f"layer {l}: no legs", flush=True); continue
    A = np.asarray(legs[l]["A"], np.float64); P = np.asarray(legs[l]["P"], np.float64)
    k = A.shape[0]
    w2b = [np.asarray(x, np.float64) for x in legs[l]["w2b"]][:k]
    eb = [np.asarray(x, np.float64) for x in legs[l]["e"]][:k]
    births = list(range(k))   # slot s was born at y_s (its P leg is W_(s+1) at layer s + 1)
    var_l = st[("var", l)]
    C_l = None
    if ("C_off", l) in st:
        C_l = st[("C_off", l)].copy(); np.fill_diagonal(C_l, var_l)
    g4 = st[("g4row", l)]
    tru = {h: mc[h]["k4"][l].astype(np.float64) for h in mc}
    r = {h: tru[h] - g4 for h in tru}
    Y = np.zeros((n, n)); Q = np.zeros((n, n)); Vt = np.zeros((n, n)); Vt1 = np.zeros((n, n))
    Ye = np.zeros((n, n)); Vte = np.zeros((n, n))   # the same with the chain's exact slice coefficient: t' = e / (2 w2)
    p4w_d = np.zeros(n); p4w_o = np.zeros(n); star_f = np.zeros(n); star_a = np.zeros(n)
    st2 = np.zeros(n); st3 = np.zeros(n); st4 = np.zeros(n); Hstk = []
    arm_err = []
    for s in range(k):
        b = births[s]
        t = st[("w1", b)] * st[("var", b)]
        At = A[s] + P[s] * t[None, :]
        X = P[s] * At * w2b[s][None, :]
        Y += X @ At.T
        Q += (At * At * w2b[s][None, :]) @ P[s].T
        Vt += X @ Lc[b]
        te = eb[s] / (2.0 * w2b[s] + 1e-30)
        Ate = A[s] + P[s] * te[None, :]
        Xe = P[s] * Ate * w2b[s][None, :]
        Ye += Xe @ At.T
        Vte += Xe @ Lc[b]
        At1 = Lc[l] @ Lc[b].T                                    # exact first-chaos cross-covariance
        arm_err.append(np.linalg.norm(At1 - At) / max(np.linalg.norm(At1), 1e-300))
        Vt1 += (P[s] * At1 * w2b[s][None, :]) @ Lc[b]
        Cb = st[("C_off", b)].copy(); np.fill_diagonal(Cb, 0.0)
        p4w_d += (X * X) @ st[("var", b)]
        p4w_o += np.sum((X @ Cb) * X, axis=1)
        mu_b, sd_b = st[("mu", b)], np.sqrt(st[("var", b)])
        al = mu_b / sd_b
        c3 = -al * phi(al) / (sd_b * sd_b)
        star_f += 4.0 * ((At ** 3) * P[s]) @ c3
        star_a += 4.0 * ((A[s] ** 3) * P[s]) @ c3
        Ps, As = P[s], A[s]
        st2 += 12.0 * ((As * As) * (Ps * Ps)) @ (c3 * t)
        st3 += 12.0 * (As * (Ps ** 3)) @ (c3 * t * t)
        st4 += 4.0 * (Ps ** 4) @ (c3 * t ** 3)
        if l in CYC_LAYERS:
            Hstk.append((b, Ps[CYC_ROWS] * w2b[s][None, :]))
    p4_ex = 12.0 * np.sum(Vt1 * Vt1, axis=1)
    p4_exc = 12.0 * np.sum(Vt * Vt, axis=1)
    p4_w = 12.0 * (p4w_d + p4w_o)
    d21c = st.get(("D21", l))
    d21m = 2.0 * Y + Q; np.fill_diagonal(d21m, 0.0)
    rel21 = np.linalg.norm(d21m - d21c) / np.linalg.norm(d21c) if d21c is not None else float("nan")
    print(f"layer {l}: k={k}  rms truth {np.sqrt(np.mean(tru['full']**2)):.3e} mean {tru['full'].mean():.3e} | chain g4 rms "
          f"{np.sqrt(np.mean(g4**2)):.3e} mean {g4.mean():.3e} | resid rms {np.sqrt(np.mean(r['full']**2)):.3e} | "
          f"D21 model rel.err {rel21:.3f} | arm rel.err (chain vs first chaos) median {np.median(arm_err):.3f} max {max(arm_err):.3f}", flush=True)
    print(f"  means: p4_exact {p4_ex.mean():.3e}  p4_exact(chain arms) {p4_exc.mean():.3e}  p4_within {p4_w.mean():.3e} "
          f"(diag part {12*p4w_d.mean():.3e})  star_full {star_f.mean():.3e}  star_V55 {star_a.mean():.3e}", flush=True)
    cand = {"p4_exact": p4_ex, "p4_exact_chainarm": p4_exc, "p4_exact_e": 12.0 * np.sum(Vte * Vte, axis=1),
            "p4_within": p4_w, "p4_within_diag": 12.0 * p4w_d,
            "p4_within_off": 12.0 * p4w_o, "star_full": star_f, "star_V55": star_a}
    cand.update({"star_A3P": star_a, "star_A2P2t": st2, "star_AP3t2": st3, "star_P4t3": st4,
                 "star_t_terms": st2 + st3 + st4})
    Ye0 = Ye.copy(); np.fill_diagonal(Ye0, 0.0); Y0 = Y.copy(); np.fill_diagonal(Y0, 0.0)
    cand["p4_dj_e"] = 12.0 * (Ye0 * Ye0) @ (1.0 / var_l)          # diagonal metric: n^2 given the hub product
    cand["p4_dj"] = 12.0 * (Y0 * Y0) @ (1.0 / var_l)
    if d21c is not None:
        cand["p4_dj_D21"] = 3.0 * (d21c * d21c) @ (1.0 / var_l)    # control: the symmetrized slice, not the physical Y
    if C_l is not None:
        Ze = np.linalg.solve(C_l + 1e-2 * var_l.mean() * np.eye(n), Ye.T)
        cand["p4_schur_e_e0.01"] = 12.0 * np.sum(Ye * Ze.T, axis=1)
        if d21c is not None:
            Zd = np.linalg.solve(C_l + 1e-2 * var_l.mean() * np.eye(n), d21c.T)
            cand["p4_schur_D21"] = 3.0 * np.sum(d21c * Zd.T, axis=1)
    if Hstk:
        # the 4-cycle 3 tr H_i^4 (first-chaos L_b) for a few neurons, against their path term
        cyc = []
        for jj, i in enumerate(CYC_ROWS):
            Hi = np.zeros((n, n))
            for (b, hrow) in Hstk:
                Hi += Lc[b].T @ (hrow[jj][:, None] * Lc[b])
            H2 = Hi @ Hi
            cyc.append((3.0 * np.sum(H2 * H2), 12.0 * float(Vt1[i] @ Vt1[i]), float(np.trace(Hi)), st[("mu", l)][i]))
        cyc = np.array(cyc)
        print(f"    4-cycle (8 neurons): 3trH^4 mean {cyc[:, 0].mean():.3e} vs path 12|HL|^2 mean {cyc[:, 1].mean():.3e}; "
              f"Euler check tr H_i vs mu_i: {np.corrcoef(cyc[:, 2], cyc[:, 3])[0, 1]:+.3f}, ratio of means {cyc[:, 2].mean() / cyc[:, 3].mean():.3f}", flush=True)
    print(f"    p4_within vs p4_exact: corr {np.corrcoef(p4_w, p4_ex)[0, 1]:+.3f}  ratio of means {p4_w.mean() / p4_ex.mean():.3f}; "
          f"chain-arm vs first-chaos arm: corr {np.corrcoef(p4_exc, p4_ex)[0, 1]:+.3f}", flush=True)
    for e in (EPS if C_l is not None else ()):
        M = C_l + e * var_l.mean() * np.eye(n)
        try:
            Z = np.linalg.solve(M, Y.T)
            cand[f"p4_schur_e{e:g}"] = 12.0 * np.sum(Y * Z.T, axis=1)
        except np.linalg.LinAlgError:
            pass
    for nm in list(cand):
        if nm.startswith("p4_schur"):
            c = cand[nm]
            print(f"    {nm}: mean {c.mean():.3e}  corr with p4_exact {np.corrcoef(c, p4_ex)[0, 1]:+.3f}  ratio of means {c.mean() / p4_ex.mean():.3f}", flush=True)
    for nm, c in cand.items():
        print(stats(nm, c, r["full"], r["h0"], r["h1"]), flush=True)
    for nm in ("p4_exact", "p4_exact_chainarm", "star_full"):
        c = cand[nm]
        print(f"    {nm:<16s} vs truth: corr(truth, c) {np.corrcoef(tru['full'], c)[0,1]:+.3f}  corr(truth, g4) {np.corrcoef(tru['full'], g4)[0,1]:+.3f}", flush=True)
    S = np.mean((r["h0"] - r["h0"].mean()) * (r["h1"] - r["h1"].mean()))
    for combo in (("p4_exact", "star_full"), ("p4_exact_e", "star_full"), ("p4_schur_e_e0.01", "star_full"),
                  ("p4_dj_e", "star_full"), ("p4_dj_e", "star_A3P", "star_t_terms"), ("p4_dj_D21", "star_full"),
                  ("p4_dj_e", "star_A3P", "star_A2P2t", "star_AP3t2", "star_P4t3")):
        if not all(c in cand for c in combo):
            continue
        Xr = np.stack([cand[c] - cand[c].mean() for c in combo], 1)
        beta, *_ = np.linalg.lstsq(Xr, r["full"] - r["full"].mean(), rcond=None)
        fit = Xr @ beta
        # split-half honesty: fit on half 0's residual, score on half 1's
        b0, *_ = np.linalg.lstsq(Xr, r["h0"] - r["h0"].mean(), rcond=None)
        r1 = r["h1"] - r["h1"].mean(); f1 = Xr @ b0
        honest = (2.0 * np.mean(r1 * f1) - np.mean(f1 * f1)) / S if S > 0 else float("nan")   # fit on h0, scored on h1
        print(f"    joint {list(combo)} slopes {' '.join(f'{x:+.3f}' for x in beta)}  signal removed {np.mean(fit * fit) / S:+.3f}"
              f"  (fit h0 -> h1: {honest:+.3f})", flush=True)
    print(f"    (residual signal rms {np.sqrt(max(S, 0)):.3e}, noise rms {np.sqrt(max(np.var(r['full']) - S, 0)):.3e})", flush=True)
np.savez(os.environ.get("C2_OUT", f"chaos2_off{net}_{tag}.npz"), out=out)
print("done", flush=True)
