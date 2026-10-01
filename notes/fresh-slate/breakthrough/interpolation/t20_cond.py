"""T20: Bolthausen-conditioned trace channel for per-neuron kappa3(z_{l+1,a}) at n=1024 (MLP 0).
Condition column w_a on its revealed m_a = mu.w_a:  w = (m/|mu|^2) mu + w_perp,  r = m/|mu|, K3 = E[(mu^.a~)^3] = t^mu.mu^:
  kappa3(z_a) ~ r^3 K3 + 3 r^2 (t^mu . w_perp) + 3 r s2 (t.mu^ - K3) + 3 s2 ((t - t^mu) . w_perp)
vs the plain chaos-1  3 s2 (t . w).  Inputs: t, t^mu, mu from t17 (exact-centred, N=131k); target kappa3 from t13 (N=2M)."""
import numpy as np
from common import bench
S = bench.load_set("w1024_d16"); W = bench.weights(S, 0).astype(np.float64); s2 = 2.0 / 1024
d = np.load("results/t17_w1024_0.npz"); h = np.load("results/t13_w1024_d16_0.npz"); N = float(h["N"]); dl = h["s1"] / N
k3 = h["s3"] / N - 3 * dl * h["s2"] / N + 2 * dl ** 3
for l in range(1, 16):
    t, tm, mu = d["t"][l - 1], d["tm"][l - 1], d["mu"][l - 1]; nm = np.linalg.norm(mu); muh = mu / nm
    Wl = W[l]; m = mu @ Wl; r = m / nm
    y = Wl.T @ t; ym = Wl.T @ tm; K3 = tm @ muh; tmu_ = t @ muh
    yperp = y - r * tmu_; ymperp = ym - r * K3
    plain = 3 * s2 * y
    cond = r ** 3 * K3 + 3 * r ** 2 * ymperp + 3 * r * s2 * (tmu_ - K3) + 3 * s2 * (yperp - ymperp)
    tgt = k3[l]; f = lambda p: 1 - ((tgt - p) ** 2).mean() / (tgt ** 2).mean()
    print(f"z layer {l+1}: explained plain {f(plain):.3f}  conditioned {f(cond):.3f}", flush=True)
