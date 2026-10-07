# The birth address of the chain's third-cumulant sources (note XXXI): how much of each read layer's D3 readout is
# carried by the birth columns of saturated birth neurons, and how a birth column's weight scales with its neuron's
# alpha at the birth layer.   python birthcols.py NET
# Source slot s is born at layer s; its legs' columns are layer-s neurons. For read layers 5, 9, 13 and thresholds
# tau, the dropped columns of slot s are those outside the top-nb by alpha_s (nb = count above tau rounded up to 64,
# exactly the V36 emulation). Reported: relative energy and correlation of the dropped part of
#   d_i = sum_s sum_j [3 w2_sj A_sij^2 P_sij + M_sij P_sij^2]
# and, per alpha bin, the mean squared column contribution sum_i t_sij^2 relative to the alpha ~ 0 bin.
import sys, os, importlib.util, numpy as np, flopscope as flops
from whestbench import MLP
net = int(sys.argv[1]); reads = (5, 9, 13)
os.environ.update({"V29_WARM_JOIN": "1", "V29_WARM_FB": "1", "V17_R_RES": "4", "V26_STRASSEN": "6", "V26_STRASSEN_MIN": "16",
                   "V32_JOIN_SMM": "1", "V32_ROT_SMM": "1", "V32_JOIN_POST": "3", "V32_JP_C": "0.1", "V21_R_OLD": "320",
                   "V24_R_OLD2": "192", "V18_R_FB": "2", "V24_AGE_OLD2": "8", "V33_K4Q": "3", "V33_K4Q_RANK": "4",
                   "V34_OPT": "abcde", "V30_DUMP_LEGS": "1", "V29_DUMP_LAYERS": ",".join(str(x) for x in range(1, 14))})
spec = importlib.util.spec_from_file_location("estv29", "est_v29.py"); mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
Wcol = np.load(f"../official/W_off{net}.npy")
mlp = MLP(width=1024, depth=16, weights=[np.ascontiguousarray(W.T).astype(np.float32) for W in Wcol])
with flops.BudgetContext(flop_budget=2**44, wall_time_limit_s=900.0, quiet=True):
    mod.Estimator().predict(mlp, 2**41)
legs = {Lg["layer"]: Lg for Lg in mod.LEGS}
alpha = {li: (Lg["mu"].astype(np.float64) / np.sqrt(Lg["var"].astype(np.float64))) for li, Lg in legs.items()}
n = 1024
def keep_cols(al, tau):
    nb = min(n, max(64, -(-int(np.sum(al > tau)) // 64) * 64))
    order = np.argsort(-al)
    m = np.zeros(n, bool); m[order[:nb]] = True
    return m
taus = (-2.5, -2.0, -1.5, -1.0, -0.5)
bins = [(-9, -2.5), (-2.5, -2.0), (-2.0, -1.5), (-1.5, -1.0), (-1.0, -0.5), (-0.5, 0.0), (0.0, 0.5), (0.5, 9)]
for li in reads:
    Lg = legs[li]
    A, P, Z, Lm = (Lg[k].astype(np.float64) for k in ("A", "P", "Z", "L"))
    w2 = np.stack(Lg["w2b"]).astype(np.float64); s = np.stack(Lg["s"]).astype(np.float64); e = np.stack(Lg["e"]).astype(np.float64)
    k = A.shape[0]
    M = P * s[:, None, :] + 3 * A * e[:, None, :] + np.einsum("kiq,kjq->kij", Z, Lm)
    t = 3 * w2[:, None, :] * A * A * P + M * P * P            # (k, n_read, n_birth)
    d = t.sum((0, 2)); E = float(d @ d)
    print(f"\n== net {net} read layer {li}: {k} sources (slot s born at layer s)", flush=True)
    for tau in taus:
        drop = np.zeros(n); fr = []
        for sl in range(1, k):
            if sl not in alpha:
                continue
            m = keep_cols(alpha[sl], tau)
            fr.append(1 - m.mean())
            drop += t[sl][:, ~m].sum(1)
        print(f"  tau {tau:+.1f}: birth columns dropped {np.mean(fr):.3f} (max {np.max(fr):.3f}); dropped part of D3: "
              f"relative energy {float(drop @ drop) / E:.2e}, corr with d {np.corrcoef(drop, d)[0, 1]:+.3f}", flush=True)
    # D21 (main legs terms): dropped part per threshold, and the readout magnitudes an oracle would need
    LA = 2 * A * P * w2[:, None, :] + P * P * e[:, None, :]
    LP = A * A * w2[:, None, :] + P * P * s[:, None, :] / 3 + (2 / 3) * M * P
    D21m = np.einsum("kij,kcj->ic", LA, A) + np.einsum("kij,kcj->ic", LP, P)
    np.fill_diagonal(D21m, 0.0); E21 = float((D21m * D21m).sum())
    for tau in taus:
        dr = np.zeros((n, n))
        for sl in range(1, k):
            if sl not in alpha:
                continue
            m = ~keep_cols(alpha[sl], tau)
            if m.any():
                dr += LA[sl][:, m] @ A[sl][:, m].T + LP[sl][:, m] @ P[sl][:, m].T
        np.fill_diagonal(dr, 0.0)
        print(f"  tau {tau:+.1f}: dropped part of D21: relative energy {float((dr * dr).sum()) / E21:.2e}", flush=True)
    var_l = Lg["var"].astype(np.float64)
    sk = Lg["D3"].astype(np.float64) / var_l ** 1.5
    print(f"  readout scales: skewness D3/var^1.5 rms {np.sqrt(np.mean(sk ** 2)):.3e}; D21 off-diag rms / var^1.5 "
          f"{np.sqrt(E21 / (n * n)) / np.mean(var_l) ** 1.5:.3e}", flush=True)
    # coefficient profiles of the birth columns against the birth alpha (all sources, pooled)
    prof = {nm: [] for nm in ("w2", "e", "s")}
    al_all = []
    for sl in range(1, k):
        if sl not in alpha:
            continue
        al_all.append(alpha[sl]); prof["w2"].append(np.abs(w2[sl])); prof["e"].append(np.abs(e[sl])); prof["s"].append(np.abs(s[sl]))
    al_all = np.concatenate(al_all)
    pr = {nm: np.concatenate(v) for nm, v in prof.items()}
    lines = []
    for lo, hi in bins:
        sel = (al_all > lo) & (al_all <= hi)
        if sel.sum() > 20:
            ref = (al_all > -0.25) & (al_all < 0.25)
            lines.append(f"[{lo:+.1f},{hi:+.1f}] " + " ".join(f"|{nm}| {pr[nm][sel].mean() / pr[nm][ref].mean():.2e}" for nm in pr))
    print("  coefficient magnitude by birth alpha (relative to |alpha| < 0.25):\n    " + "\n    ".join(lines), flush=True)
    # column weight against the birth neuron's alpha
    rows = []
    for sl in range(1, k):
        if sl not in alpha:
            continue
        cw = (t[sl] ** 2).sum(0)                          # (n_birth,)
        rows.append((alpha[sl], cw / np.median(cw)))
    al = np.concatenate([r[0] for r in rows]); cw = np.concatenate([r[1] for r in rows])
    ref = cw[(al > -0.25) & (al < 0.25)].mean()
    phi = lambda a: np.exp(-a * a / 2) / np.sqrt(2 * np.pi)
    out = []
    for lo, hi in bins:
        sel = (al > lo) & (al <= hi)
        if sel.sum() > 20:
            amid = al[sel].mean()
            out.append(f"[{lo:+.1f},{hi:+.1f}] n={sel.sum()} w={cw[sel].mean() / ref:.2e} (phi^2 ratio {(phi(amid) / phi(0)) ** 2:.2e})")
    print("  column weight sum_i t_ij^2 by birth alpha (relative to |alpha| < 0.25):\n    " + "\n    ".join(out), flush=True)
