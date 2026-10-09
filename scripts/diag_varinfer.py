"""Infer the per-neuron pre-activation variance that the readout formula would need, given the truth means of the
previous layer (exact to 1e-5) and the chain's own kappa3/kappa4, and compare with the chain's variance. Then regress
the relative variance error on the loadings of the top eigenmodes of the chain's post-activation covariance.
  python scripts/diag_varinfer.py DATA NET [opt=val ...]"""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from scipy.special import ndtr
from whest.k3chain2 import k3_chain2
from whest.relu_gauss import phi
D, net = sys.argv[1], int(sys.argv[2]); opts = {}
for a in sys.argv[3:]:
    k, v = a.split("="); opts[k] = v if k == "k4" else (float(v) if "." in v else int(v))
W = np.load(f"{D}/W_off{net}.npy").astype(np.float64); mt = np.load(f"{D}/truth_off{net}.npz")["m"].astype(np.float64)
rec = {}; out, fl = k3_chain2(W, opts, record=rec); L, n, _ = W.shape
def readout(mu, var, D3, K4):
    s = np.sqrt(var); a = mu / s; p = phi(a)
    return s * (p + a * ndtr(a)) - D3 * a * p / (6 * var) + K4 * (a * a - 1) * p / (24 * s ** 3)
print("layer | rms rel var err (inferred vs chain) | mean (coherent) | R^2 on top-1/4/16 mode loadings | corr(err, loading1) | |resid| after top-16")
for l in range(1, L):
    r = rec[l]; mu = W[l] @ mt[l - 1]; D3 = r["D3"]; K4 = r.get("K4", np.zeros(n)) if "K4" in r else np.zeros(n)
    v0 = r["var"].copy(); lo = v0 * 0.5; hi = v0 * 2.0
    for _ in range(60):   # bisection in var (readout increasing in var for the Gaussian part)
        mid = np.sqrt(lo * hi); f = readout(mu, mid, D3, K4) - mt[l]
        lo = np.where(f < 0, mid, lo); hi = np.where(f >= 0, mid, hi)
    vinf = np.sqrt(lo * hi); err = (v0 - vinf) / vinf
    # loadings: eigenmodes of the chain's previous post-activation covariance Kh_{l-1} = rec[l]["C"] is pre-activation C_l;
    # use C_l's own eigenmodes: fraction of var_i in mode k = lambda_k q_ki^2 / var_i
    C = r["C"]; w, Q = np.linalg.eigh(C); idx = np.argsort(w)[::-1]; w = w[idx]; Q = Q[:, idx]
    load = (Q ** 2) * w[None, :] / v0[:, None]      # (n, modes)
    def r2(k):
        X = np.c_[np.ones(n), load[:, :k]]; beta, *_ = np.linalg.lstsq(X, err, rcond=None); res = err - X @ beta
        return 1 - np.sum(res ** 2) / np.sum((err - err.mean()) ** 2), np.sqrt(np.mean(res ** 2))
    print(f"{l:5d} | {np.sqrt(np.mean(err**2)):8.2e} | {err.mean():+8.2e} | {r2(1)[0]:5.3f} {r2(4)[0]:5.3f} {r2(16)[0]:5.3f} | {np.corrcoef(err, load[:,0])[0,1]:+6.3f} | {r2(16)[1]:8.2e}   top eigenvalue share of tr C: {w[0]/w.sum():.3f}")
