# Kernel-level test of the conditional-contraction development (DEVELOPMENT.md) on the chain's real legs.
#   python cohere.py NET LAYERS(comma list)
# The chain's D3 readout is a cubic contraction over the internal index a = (source k, birth neuron j):
#   D3_i ~ sum_a X1_ia X2_ia X3_ia,  t1 = (A, A, 3 w2 P),  t2 = (M, P, P),  M = P s + 3 A e + Z L^T
# (the feedback and feed terms are left out and reported as a share). For it we compute, per output row i:
#   kappa_i = |d_i| / sum_a |t_ia|                      (sign coherence; no positive gauge changes t_ia)
#   raw cubic, ungauged (Thm 2): ||X1_i||^2 ||X2_i||^2 ||X3_i||^2 + 3 d_i^2 - 4 sum t^2
#   raw cubic, best per-row gauge (Holder (8)): (sum |t|^{2/3})^3 + 3 d^2 - 4 sum t^2
#   each filling (A|BC, B|AC, C|AB), ungauged (Thm 3 at u = 0), with the convex SHARED gauge optimum, and the per-row
#   floor (sum |t|)^2 + d^2 - 2 sum t^2 (identical for all three cuts); shared / per-row = the cycle (Hodge) loss.
# All variances are single-probe, summed over rows, and reported relative to sum_i d_i^2 (so m probes give relative
# squared error ratio/m). D21 (matrix readout, sum_a LA_ia A_ca + LP_ia P_ca) gets kappa and the per-entry floors on
# a sample of entries.
import sys, os, importlib.util, numpy as np, flopscope as flops
from scipy.optimize import minimize
from whestbench import MLP
net = int(sys.argv[1]); layers = sys.argv[2]
os.environ.update({"V29_WARM_JOIN": "1", "V29_WARM_FB": "1", "V17_R_RES": "4", "V26_STRASSEN": "6", "V26_STRASSEN_MIN": "16",
                   "V32_JOIN_SMM": "1", "V32_ROT_SMM": "1", "V32_JOIN_POST": "3", "V32_JP_C": "0.1", "V21_R_OLD": "320",
                   "V24_R_OLD2": "192", "V18_R_FB": "2", "V24_AGE_OLD2": "8", "V33_K4Q": "3", "V33_K4Q_RANK": "4",
                   "V34_OPT": "abcde", "V30_DUMP_LEGS": "1", "V29_DUMP_LAYERS": layers})
spec = importlib.util.spec_from_file_location("estv29", "est_v29.py"); mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
Wcol = np.load(f"../official/W_off{net}.npy")
mlp = MLP(width=1024, depth=16, weights=[np.ascontiguousarray(W.T).astype(np.float32) for W in Wcol])
with flops.BudgetContext(flop_budget=2**44, wall_time_limit_s=900.0, quiet=True):
    mod.Estimator().predict(mlp, 2**41)


def gauge_opt(X2, Y2, maxiter=400):
    # min_u sum_i R_i S_i, R = X2 e^{2u}, S = Y2 e^{-2u}  (convex; matrix-free gradient (12)); bounded range |u| <= 8
    r = X2.shape[1]
    def f(u):
        ep, em = np.exp(2 * u), np.exp(-2 * u)
        R, S = X2 @ ep, Y2 @ em
        val = float(R @ S)
        g = 2 * (ep * (X2.T @ S) - em * (Y2.T @ R))
        return val, g
    res = minimize(f, np.zeros(r), jac=True, method="L-BFGS-B", bounds=[(-8, 8)] * r, options=dict(maxiter=maxiter))
    return res.fun, res.x, res.nit


rng = np.random.default_rng(0)
for Lg in mod.LEGS:
    li = Lg["layer"]
    A, P, Z, Lm = (Lg[k].astype(np.float64) for k in ("A", "P", "Z", "L"))
    w2 = np.stack(Lg["w2b"]).astype(np.float64); s = np.stack(Lg["s"]).astype(np.float64); e = np.stack(Lg["e"]).astype(np.float64)
    k, n, _ = A.shape
    M = P * s[:, None, :] + 3 * A * e[:, None, :] + np.einsum("kiq,kjq->kij", Z, Lm)
    # internal index a = (term, k, j): move it last
    X1 = np.concatenate([A, M], 0).transpose(1, 0, 2).reshape(n, -1)
    X2 = np.concatenate([A, P], 0).transpose(1, 0, 2).reshape(n, -1)
    X3 = np.concatenate([3 * w2[:, None, :] * P, P], 0).transpose(1, 0, 2).reshape(n, -1)
    t = X1 * X2 * X3
    d = t.sum(1); D3 = Lg["D3"].astype(np.float64)
    at = np.abs(t).sum(1); t2s = (t * t).sum(1); sd2 = (d * d).sum()
    kap = np.abs(d) / at
    raw_ung = ((X1 * X1).sum(1) * (X2 * X2).sum(1) * (X3 * X3).sum(1) + 3 * d * d - 4 * t2s).sum() / sd2
    raw_hold = ((np.abs(t) ** (2 / 3)).sum(1) ** 3 + 3 * d * d - 4 * t2s).sum() / sd2
    floor = (at * at + d * d - 2 * t2s).sum() / sd2
    print(f"\n== net {net} layer {li}: k={k} sources (ka={Lg['ka']} old), internal rank r={X1.shape[1]}; "
          f"main+M terms carry corr {np.corrcoef(d, D3)[0, 1]:.3f}, energy share {np.sum(d * d) / np.sum(D3 * D3):.3f} of D3", flush=True)
    print(f"  sign coherence kappa=|d|/sum|t|: median {np.median(kap):.4f}  mean {kap.mean():.4f}  "
          f"(random signs over r terms ~ {np.sqrt(2 / np.pi / X1.shape[1]):.4f})")
    print(f"  single-probe variance / sum d^2:  raw ungauged {raw_ung:.3e} | raw best per-row gauge (Holder) {raw_hold:.3e}"
          f" | filled per-row floor {floor:.3e}")
    for name, X, Y in (() if os.environ.get("COH_SKIPG") == "1" else
                       (("A|BC", X1, X2 * X3), ("B|AC", X2, X1 * X3), ("C|AB", X3, X1 * X2))):
        Xs, Ys = X * X, Y * Y
        ung = ((Xs.sum(1) * Ys.sum(1)) + d * d - 2 * t2s).sum() / sd2
        G, u, nit = gauge_opt(Xs, Ys)
        shared = (G + (d * d - 2 * t2s).sum()) / sd2
        print(f"  {name}: ungauged {ung:.3e} | shared gauge {shared:.3e} (L-BFGS {nit} it, |u| max {np.abs(u).max():.2f})"
              f" | cycle loss shared/floor {shared / floor:.3f}", flush=True)
    # where the cancellation lives: within a source (across birth neurons j) or across sources / terms
    tk = t.reshape(n, 2 * k, n)                     # (row, term-source, j)
    dk = tk.sum(2); ak = np.abs(tk).sum(2)          # per term-source net and absolute mass
    kin = np.abs(dk) / np.maximum(ak, 1e-300)
    kacross = np.abs(dk.sum(1)) / np.abs(dk).sum(1)
    share_new = np.abs(dk[:, [k - 1, 2 * k - 1]]).sum(1) / np.abs(dk).sum(1)
    srt = -np.sort(-np.abs(t), axis=1); cum = np.cumsum(srt, 1) / at[:, None]
    hq = {h: float(np.median(cum[:, h - 1])) for h in (16, 128, 1024)}
    print(f"  cancellation: within a term-source over j kappa median {np.median(kin):.4f}; across the {2 * k} term-sources "
          f"kappa median {np.median(kacross):.3f}; newest source's share of sum_k |d_k| median {np.median(share_new):.3f}; "
          f"abs-mass share of the top 16/128/1024 terms per row {hq[16]:.3f}/{hq[128]:.3f}/{hq[1024]:.3f}", flush=True)
    # D21 (matrix readout): per-entry coherence on a sample of (i, c) with i != c
    w2c = w2[:, None, :]; ec = e[:, None, :]; sc = s[:, None, :]
    LA = 2 * A * P * w2c + P * P * ec
    LP = A * A * w2c + P * P * sc / 3 + (2 / 3) * M * P
    D21 = Lg["D21"].astype(np.float64)
    ii = rng.integers(0, n, 4096); cc = rng.integers(0, n, 4096); keep = ii != cc; ii, cc = ii[keep], cc[keep]
    T21 = np.concatenate([LA[:, ii, :] * A[:, cc, :], LP[:, ii, :] * P[:, cc, :]], 0)   # (2k, m, n): terms of each entry
    d21 = T21.sum((0, 2)); a21 = np.abs(T21).sum((0, 2)); q21 = (T21 * T21).sum((0, 2))
    kap21 = np.abs(d21) / a21
    fl21 = (a21 * a21 + d21 * d21 - 2 * q21).sum() / (d21 * d21).sum()
    print(f"  D21 hub terms: corr with dumped D21 {np.corrcoef(d21, D21[ii, cc])[0, 1]:.3f}; kappa median {np.median(kap21):.4f}; "
          f"filled per-entry floor / sum d^2 {fl21:.3e} (probes for 10% rel. error >= {100 * fl21:.2e})", flush=True)
