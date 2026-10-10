"""Aggregate exact-leaf results: python scripts/leaf_stats.py DATA NET FILE [FILE ...]
Reports events per leaf (vs 2nL), the wall-identity residuals, the plane-mean estimate c_n mean_P g against the
truth, and the per-neuron transverse variance c_n^2 Var_P(g_i) against the neuron variance (avg_variance)."""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from whest.leaf import c_radial

D, net = sys.argv[1], int(sys.argv[2]); files = sys.argv[3:]
t = np.load(f"{D}/truth_off{net}.npz"); m = t["m"][-1].astype(np.float64); v_avg = float(t["avg_variance"])
G = np.concatenate([np.load(f)["g"] for f in files]); E = np.concatenate([np.load(f)["events"] for f in files])
R = np.concatenate([np.load(f)["wall"] for f in files]); S = np.concatenate([np.load(f)["secs"] for f in files])
K, n = G.shape; L = E.shape[1]; cn = c_radial(n)
Y = cn * G                                                   # per-plane unbiased estimates of the means
print(f"net {net}: {K} leaves; events per leaf {E.sum(1).mean():.0f} (2nL = {2*n*L}); per layer mean {np.round(E.mean(0)).astype(int).tolist()}")
print(f"wall residual max {R.max():.1e}; seconds per leaf {S.mean():.0f}")
est = Y.mean(0); se2 = Y.var(0, ddof=1) / K
print(f"plane-mean estimate vs truth: MSE {np.mean((est - m)**2):.3e} (predicted from the plane variance {se2.mean():.3e}); truth rms {np.sqrt(np.mean(m**2)):.4f}")
tv = Y.var(0, ddof=1)
print(f"per-neuron transverse variance c_n^2 Var_P(g): mean {tv.mean():.4f}, median {np.median(tv):.4f}; neuron variance avg_variance = {v_avg:.4f}; ratio mean {tv.mean()/v_avg:.3f}")
print(f"  (Theorem 5: ratio = sum_r lambda_2r ||f_2r||^2 / Var f with lambda_2 = {(n-2)/(2*(n-1)):.3f}; a ratio near 0.5 means degree 2 dominates, near 0 means high degrees)")
print(f"per-leaf estimate rms error {np.sqrt(np.mean((Y - m)**2)):.4f}; leaves needed for MSE 1e-8 at this variance: {tv.mean()/1e-8:.1e}")
