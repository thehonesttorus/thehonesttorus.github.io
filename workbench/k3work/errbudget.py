# The chain's final error as a linear functional of its per-layer statistic errors (note XXXI).
#   python errbudget.py NET MCPREFIX        (needs chain_off{NET}.npz from oracle_one.py NET dump)
# The mean map of layer l reads only the marginal law of each z_li: E relu(z) = mu Phi(a) + sigma phi(a) +
# sigma sum_k kappa_k / (k! sigma^k) He_(k-2)(-a) phi(a) (Gram-Charlier), so to first order an error of the variance,
# third or fourth cumulant diagonal injects into E[y_l]
#   var: phi / (2 sigma) d_var      kappa_3: -a phi / (6 sigma^2) d_k3      kappa_4: (a^2 - 1) phi / (24 sigma^3) d_k4
# and an error of E[y_(l-1)] reaches E[y_l] through Phi(a_l) o (W_l .). Everything else the chain carries (C_off, D21, the
# kappa_4 slices, the sources) acts on the output only through these three diagonals at later layers. With the
# Jacobian evaluated at Monte Carlo truth and the injections s_l^X from the chain's dumped statistics minus truth,
# Delta^X = sum_l J_(15 <- l) s_l^X is the first-order part of the output error carried by channel X, and the oracle that
# replaces X at every layer is predicted to leave e - Delta^X (e = the chain's output error). Noise-free estimates use
# the two Monte Carlo halves: <e, Delta> from their mean, |Delta|^2 from <Delta_h0, Delta_h1>.
import sys, os, math, importlib.util, numpy as np, flopscope as flops
from whestbench import MLP
net, pre = int(sys.argv[1]), sys.argv[2]
os.environ.update({"V29_WARM_JOIN": "1", "V29_WARM_FB": "1", "V17_R_RES": "4", "V26_STRASSEN": "6", "V26_STRASSEN_MIN": "16",
                   "V32_JOIN_SMM": "1", "V32_ROT_SMM": "1", "V32_JOIN_POST": "3", "V32_JP_C": "0.1", "V21_R_OLD": "320",
                   "V24_R_OLD2": "192", "V18_R_FB": "2", "V24_AGE_OLD2": "8", "V33_K4Q": "3", "V33_K4Q_RANK": "4",
                   "V34_OPT": "abcde"})
spec = importlib.util.spec_from_file_location("estv29", "est_v29.py"); mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
Wcol = np.load(f"../official/W_off{net}.npy").astype(np.float64); mt = np.load(f"../official/truth_off{net}.npz")["m"].astype(np.float64)
mlp = MLP(width=1024, depth=16, weights=[np.ascontiguousarray(W.T).astype(np.float32) for W in Wcol])
with flops.BudgetContext(flop_budget=2**44, wall_time_limit_s=900.0, quiet=True):
    out = np.asarray(mod.Estimator().predict(mlp, 2**41), dtype=np.float64)
e = out[-1] - mt[-1]; E = float(e @ e); n = e.size
ch = np.load(os.environ.get("CHAIN_DUMP", f"chain_off{net}.npz")); M = {t: np.load(f"{pre}_{t}.npz") for t in ("full", "h0", "h1")}
F = M["full"]; L = Wcol.shape[0]
mu, var = F["mu"].astype(np.float64), F["var"].astype(np.float64)
sig = np.sqrt(var); a = mu / sig; ph = np.exp(-a * a / 2) / math.sqrt(2 * math.pi); Ph = 0.5 * (1 + np.vectorize(math.erf)(a / math.sqrt(2)))
coef = {"var": ph / (2 * sig), "k3": -a * ph / (6 * var), "k4": (a * a - 1) * ph / (24 * var * sig)}
cname = {"var": "var", "k3": "D3", "k4": "g4row"}


def prop(l, v):
    for k in range(l + 1, L):
        v = Ph[k] * (Wcol[k] @ v)
    return v


D = {}          # D[(X, truth)] = propagated vector; P[(X, truth)] = per-layer propagated vectors
P = {}
for X in ("var", "k3", "k4"):
    for t in ("full", "h0", "h1"):
        tot = np.zeros(n); per = {}
        for l in range(1, L):
            key = f"{cname[X]}_{l}"
            if key not in ch.files:
                continue
            s = coef[X][l] * (ch[key].astype(np.float64) - M[t][X][l].astype(np.float64))
            per[l] = prop(l, s); tot += per[l]
        D[(X, t)] = tot; P[(X, t)] = per
print(f"net {net}: chain output MSE {E / n:.4e} (vs official truth); first-order channel content, noise-free:")
combos = [("var",), ("k3",), ("k4",), ("k3", "k4"), ("var", "k3", "k4")]
for c in combos:
    d0 = sum(D[(X, "h0")] for X in c); d1 = sum(D[(X, "h1")] for X in c); df = sum(D[(X, "full")] for X in c)
    cross = 0.5 * float(e @ d0 + e @ d1); size = float(d0 @ d1)
    pred = (E - 2 * cross + size) / E - 1
    print(f"  {'+'.join(c):10s}: predicted oracle change {100 * pred:+6.1f}% | |Delta|^2/|e|^2 {size / E:.3f} "
          f"(with MC noise at full N {float(df @ df) / E:.3f}) | corr(e, Delta) {cross / np.sqrt(E * max(size, 1e-300)):+.3f}")
d0 = sum(D[(X, "h0")] for X in ("var", "k3", "k4")); d1 = sum(D[(X, "h1")] for X in ("var", "k3", "k4"))
r0, r1 = e - d0, e - d1
print(f"  unexplained by the three channels to first order: |e - Delta|^2/|e|^2 = {float(r0 @ r1) / E:.3f}")
print("  per source layer, noise-free |J s_l|^2 / |e|^2 and its correlation with e (var | k3 | k4):")
for l in range(1, L):
    row = []
    for X in ("var", "k3", "k4"):
        if l in P[(X, "h0")]:
            p0, p1 = P[(X, "h0")][l], P[(X, "h1")][l]
            sz = float(p0 @ p1)
            row.append(f"{sz / E:6.3f} ({0.5 * float(e @ (p0 + p1)) / np.sqrt(E * max(abs(sz), 1e-300)):+.2f})")
    print(f"    layer {l:2d}: " + " | ".join(row), flush=True)
