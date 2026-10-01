"""T22: is the bulk conditionally Gaussian given the collective coordinate g_l = e1 . (a_l - mu) (e1 = top eigvector of
Cov(a_l), from t15)?  MC at n=1024 (MLP 0): per-neuron kappa3, kappa4 of z_{l+1} unconditional vs pooled within-bin
(8 quantile bins of g_l); also kappa3(g), kappa4(g) and the share of per-neuron kappa3 explained by the trace channel."""
import numpy as np
from common import bench
S = bench.load_set("w1024_d16"); W = bench.weights(S, 0); T = S["means"][0]; n = 1024
C = np.load("/tmp/claude-0/s/t15_w1024_0.npz")["C"]
layers = [4, 9, 13]; e1 = {l: np.linalg.eigh(C[l].astype(np.float64))[1][:, -1].astype(np.float32) for l in layers}
for l in layers:
    if e1[l] @ T[l] < 0: e1[l] = -e1[l]
    print(f"layer {l+1}: top-eigvec share of Cov {np.linalg.eigvalsh(C[l].astype(np.float64))[-1]/np.trace(C[l]):.3f}, cos(e1, mu) {e1[l]@T[l]/np.linalg.norm(T[l]):.3f}")
N = 262144; chunk = 8192; rng = np.random.default_rng(77)
G = {l: [] for l in layers}; Z = {l: [] for l in layers}
for c in range(N // chunk):
    h = rng.standard_normal((chunk, n)).astype(np.float32)
    for l in range(max(layers) + 2):
        z = h @ W[l]
        if l - 1 in layers: Z[l - 1].append(z[:, :256].copy())     # 256 neurons of z_{l+1}
        h = np.maximum(z, 0.0)
        if l in layers: G[l].append((h - T[l].astype(np.float32)) @ e1[l])
for l in layers:
    g = np.concatenate(G[l]).astype(np.float64); z = np.concatenate(Z[l]).astype(np.float64)
    gs = (g - g.mean()) / g.std()
    k3g = (gs ** 3).mean(); k4g = (gs ** 4).mean() - 3
    y = z - z.mean(0); v = (y ** 2).mean(0)
    k3u = (y ** 3).mean(0) / v ** 1.5; k4u = (y ** 4).mean(0) / v ** 2 - 3
    bins = np.quantile(g, np.linspace(0, 1, 9)); idx = np.clip(np.searchsorted(bins, g) - 1, 0, 7)
    k3c = np.zeros(256); k4c = np.zeros(256)
    for b in range(8):
        yb = z[idx == b]; yb = yb - yb.mean(0); vb = (yb ** 2).mean(0)
        k3c += (yb ** 3).mean(0) / vb ** 1.5 / 8; k4c += ((yb ** 4).mean(0) / vb ** 2 - 3) / 8
    # binning a Gaussian into 8 quantile bins itself changes within-bin shape only along g; report rms
    r = lambda x: np.sqrt((x ** 2).mean())
    print(f"z layer {l+2}: kappa3(g) {k3g:.2f} kappa4(g) {k4g:.2f} | per-neuron std. skew rms uncond {r(k3u):.4f} within-g {r(k3c):.4f} | "
          f"excess kurt mean uncond {k4u.mean():+.4f} within-g {k4c.mean():+.4f} (rms {r(k4u):.4f} / {r(k4c):.4f}) | corr(z_j, g) rms {r((y*gs[:,None]).mean(0)/np.sqrt(v)):.3f}", flush=True)
