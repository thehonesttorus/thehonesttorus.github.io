"""T13: high-N single-pass MC (centred on bake means; bake noise << signal) of per-neuron z variance and third
central moment at every layer, for the oracle substitutions of t12.  Saves results/t13_<name>_<i>_N.npz"""
import sys, numpy as np, time
from common import bench
name, i, N = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]); chunk = 8192
S = bench.load_set(name); W = bench.weights(S, i); T = S["means"][i]; L, n, _ = W.shape
mz = np.stack([np.zeros(n)] + [T[l - 1] @ W[l].astype(np.float64) for l in range(1, L)]).astype(np.float32)
s1 = np.zeros((L, n)); s2 = np.zeros((L, n)); s3 = np.zeros((L, n)); rng = np.random.default_rng(2024); done = 0; t0 = time.time()
while done < N:
    h = rng.standard_normal((chunk, n)).astype(np.float32)
    for l in range(L):
        z = h @ W[l]; y = z - mz[l]
        s1[l] += y.sum(0, dtype=np.float64); y2 = y * y; s2[l] += y2.sum(0, dtype=np.float64); s3[l] += (y2 * y).sum(0, dtype=np.float64)
        h = np.maximum(z, 0.0)
    done += chunk
    if done % (chunk * 32) == 0:
        print(done, f"{time.time()-t0:.0f}s", flush=True)
        np.savez(f"results/t13_{name}_{i}.npz", s1=s1, s2=s2, s3=s3, N=done, mz=mz)
np.savez(f"results/t13_{name}_{i}.npz", s1=s1, s2=s2, s3=s3, N=done, mz=mz); print("done", flush=True)
