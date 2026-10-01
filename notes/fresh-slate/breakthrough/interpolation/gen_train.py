"""Training networks for disorder-conditional corrections: fresh seeds (not bench seeds), MC truth of all-layer means.
Usage: gen_train.py width depth seed0 count N  -> results/train/w{width}_s{seed}.npz (means float64 (L,n), N, var_final)"""
import sys, os, time
import numpy as np
from common import bench
w, d, s0, cnt, N = map(int, sys.argv[1:6])
os.makedirs("results/train", exist_ok=True)
for s in range(s0, s0 + cnt):
    f = f"results/train/w{w}_d{d}_s{s}.npz"
    if os.path.exists(f): continue
    t0 = time.time()
    W = bench.weights_from_seed(s, w, d)
    rng = np.random.default_rng(10_000 + s)
    acc = np.zeros((d, w)); acc2 = np.zeros(w); chunk = 8192; done = 0
    while done < N:
        h = rng.standard_normal((chunk, w)).astype(np.float32)
        for l in range(d):
            h = np.maximum(h @ W[l], 0.0); acc[l] += h.sum(0, dtype=np.float64)
        acc2 += (h.astype(np.float64) ** 2).sum(0); done += chunk
    mean = acc / done; var = acc2 / done - mean[-1] ** 2
    np.savez(f, means=mean, N=done, noise=var.mean() / done)
    print(s, f"{time.time()-t0:.0f}s noise {var.mean()/done:.2e}", flush=True)
