"""T14: is the non-trace (chaos>=2) residual of per-neuron kappa3(z_{l+1}) the age-0 Gaussian source?
residual r = kappa3_MC(z_{l+1}) - 3 sigma^2 W^T t_MC(l);  age-0 source g = kappa3 of (relu(Z) W_{l+1}), Z ~ N(m_l, S_l)
(TC state), minus its own chaos-1 part 3 sigma^2 W^T t_G with t_G from the same Gaussian samples."""
import glob, numpy as np
from common import bench
import tc
S_ = bench.load_set("w1024_d16"); W = bench.weights(S_, 0); n = 1024; sig2 = 2.0 / n
d = dict(np.load(sorted(glob.glob("results/t9_w1024_d16_0_N*.npz"))[-1]))
import os
if os.path.exists("results/t13_w1024_d16_0.npz"):
    h = np.load("results/t13_w1024_d16_0.npz"); Nh = float(h["N"]); dl = h["s1"] / Nh
    d["m3"] = h["s3"] / Nh - 3 * dl * h["s2"] / Nh + 2 * dl ** 3; print("using t13 kappa3, N =", Nh)
_, states, ts = tc.predict(W, ret_C="states")
rng = np.random.default_rng(5); N = 131072; chunk = 8192
for l in (3, 7, 11, 13):
    m, S = states[l]; Lc = np.linalg.cholesky(S + 1e-6 * np.trace(S) / n * np.eye(n)).astype(np.float32)
    A = []
    for c in range(N // chunk):
        Z = rng.standard_normal((chunk, n)).astype(np.float32) @ Lc.T + m.astype(np.float32)
        A.append(np.maximum(Z, 0))
    A = np.concatenate(A).astype(np.float64); A -= A.mean(0)
    tG = ((A ** 2).sum(1)[:, None] * A).mean(0)
    Y = A @ W[l + 1].astype(np.float64); k3G = (Y ** 3).mean(0); del A, Y
    g = k3G - 3 * sig2 * (W[l + 1].astype(np.float64).T @ tG)
    r = d["m3"][l + 1] - 3 * sig2 * (W[l + 1].astype(np.float64).T @ d["t"][l])
    k3 = d["m3"][l + 1]
    print(f"z layer {l+2}: rms k3 {np.sqrt((k3**2).mean()):.4f} rms resid {np.sqrt((r**2).mean()):.4f} rms age0-chaos3 {np.sqrt((g**2).mean()):.4f} "
          f"corr(resid, age0) {np.corrcoef(r, g)[0,1]:.3f}  frac of resid explained {1-((r-g)**2).mean()/(r**2).mean():.3f}", flush=True)
