"""Pool the circ jobs: python scripts/circle_stats.py RESULTS_DIR DATA NET
Per layer: per-neuron variance (averaged over neurons) of a single Gaussian sample, an antithetic pair, an exact
great-circle leaf (c_n g) and a rotated cross-polytope design; each estimator's mean against the truth (z-score of
the bias); and the harmonic split implied by Theorem 5 of the circle note:
  V_anti = e2 + e4 + ...,   V_leaf = lam2 e2 + lam4 e4 + ...   =>  e2 ~ (V_leaf - lam4 V_anti)/(lam2 - lam4)."""
import sys, os, glob, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from whest.leaf import c_radial

R, D, net = sys.argv[1], sys.argv[2], int(sys.argv[3])
truth = np.load(f"{D}/truth_off{net}.npz")["m"].astype(np.float64); L, n = truth.shape; cn = c_radial(n)
al = (n - 2) / 2; lam2 = (n - 2) / (2 * (n - 1)); lam4 = lam2 * (al + 1) * 1.5 / (2 * (al + 1.5))
pool = {}
for f in glob.glob(f"{R}/harm_{net}_*.npz"):
    z = np.load(f)
    for k in ("single", "antithetic", "crosspoly"):
        p = pool.setdefault(k, [0.0, 0.0, 0]); p[0] = p[0] + z[f"{k}_0"]; p[1] = p[1] + z[f"{k}_1"]; p[2] += int(z[f"{k}_2"])
G = np.concatenate([np.load(f)["g"] for f in sorted(glob.glob(f"{R}/leafb_{net}_*.npz"))]) * cn      # (K, L, n)
secs = np.concatenate([np.load(f)["secs"] for f in glob.glob(f"{R}/leafb_{net}_*.npz")])
walls = np.concatenate([np.load(f)["walls"] for f in glob.glob(f"{R}/leafb_{net}_*.npz")])
K = G.shape[0]
print(f"net {net}: {pool['single'][2]} pairs, {pool['crosspoly'][2]} designs, {K} leaves (walls/leaf {walls.sum(1).mean():.0f}, {secs.mean():.1f}s/leaf); lam2 {lam2:.4f} lam4 {lam4:.4f}")
print(" l |  V_single  V_anti  V_leaf  V_cross  | bias z: single anti leaf cross | odd share  e2 share  e>=4 share (of V_single)")
for l in range(L):
    V = {}; Mn = {}
    for k, (s1, s2, c) in pool.items():
        mu = s1[l] / c; V[k] = s2[l] / c - mu * mu; Mn[k] = (mu, c)
    vl = G[:, l, :].var(0, ddof=1); ml = G[:, l, :].mean(0)
    z = lambda mu, var, c: np.mean((mu - truth[l]) ** 2) / np.mean(var / c)       # ~1 if unbiased
    vs, va, vc = V["single"].mean(), V["antithetic"].mean(), V["crosspoly"].mean()
    e2 = (vl.mean() - lam4 * va) / (lam2 - lam4); e4 = va - e2
    print(f"{l:2d} | {vs:.2e} {va:.2e} {vl.mean():.2e} {vc:.2e} | {z(Mn['single'][0], V['single'], Mn['single'][1]):.2f} {z(Mn['antithetic'][0], V['antithetic'], Mn['antithetic'][1]):.2f} "
          f"{z(ml, vl, K):.2f} {z(Mn['crosspoly'][0], V['crosspoly'], Mn['crosspoly'][1]):.2f} | {1 - 2 * va / vs if False else (vs - 2 * va) / vs:.3f}  {e2 / vs:.3f}  {e4 / vs:.3f}")
print("(odd share = (V_single - 2 V_anti)/V_single, since f(x) and f(-x) share the even part: V_single = odd + even, V_anti = even.)")
