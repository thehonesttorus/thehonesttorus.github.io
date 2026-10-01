"""T9: exact-centring MC of t_l = E[|a-mu|^2 (a-mu)] (two passes over identical samples) and per-neuron kappa3(z),
at n=1024; compares with the TC self-consistent t and the chaos-1 prediction 3 sigma^2 W^T t."""
import sys, numpy as np
from common import bench
import tc
name, i, N = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]); chunk = 8192
S = bench.load_set(name); W = bench.weights(S, i); L, n, _ = W.shape; sig2 = 2.0 / n
def passes(fn):
    rng = np.random.default_rng(99)
    for c in range(N // chunk):
        h = rng.standard_normal((chunk, n)).astype(np.float32)
        for l in range(L):
            z = h @ W[l]; h = np.maximum(z, 0.0); fn(l, z, h)
sa = np.zeros((L, n)); sz = np.zeros((L, n))
def p1(l, z, h): sa[l] += h.sum(0, dtype=np.float64); sz[l] += z.sum(0, dtype=np.float64)
passes(p1); mu = sa / N; mz = sz / N
t = np.zeros((L, n)); m2 = np.zeros((L, n)); m3 = np.zeros((L, n))
def p2(l, z, h):
    u = h.astype(np.float64) - mu[l]; t[l] += ((u ** 2).sum(1)[:, None] * u).sum(0)
    y = z.astype(np.float64) - mz[l]; m2[l] += (y ** 2).sum(0); m3[l] += (y ** 3).sum(0)
passes(p2); t /= N; m2 /= N; m3 /= N
np.savez(f"results/t9_{name}_{i}_N{N}.npz", t=t, mu=mu, mz=mz, m2=m2, m3=m3)
p, ts = tc.predict(W, ret_state=True)
for l in range(L):
    r = ts[l] @ t[l] / (t[l] @ t[l]); c = np.corrcoef(ts[l], t[l])[0, 1]
    msg = f"layer {l+1}: t_tc vs t_mc corr {c:.4f} proj ratio {r:.3f}"
    if l > 0:
        c1 = 3 * sig2 * (t[l - 1] @ W[l].astype(np.float64)); k3 = m3[l]
        msg += f" | kappa3 chaos-1 frac {1 - ((k3 - c1) ** 2).mean() / (k3 ** 2).mean():.3f}"
    print(msg, flush=True)
