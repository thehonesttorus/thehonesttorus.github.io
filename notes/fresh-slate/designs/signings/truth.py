"""Monte Carlo truth for Stage Q (weights from whestbench.generation.sample_mlp(seed)).
Accumulates per-layer sums and final-layer sum of squares in float64; resumable in chunks.
usage: python truth.py WIDTH SEED NSAMPLES OUTDIR"""
import sys, os, time, numpy as np
from whestbench.generation import sample_mlp

def weights(width, depth, seed):
    m = sample_mlp(width, depth, rng=np.random.default_rng(seed), seed=seed) if 'rng' in sample_mlp.__code__.co_varnames else sample_mlp(width, depth, seed=seed)
    return np.stack([np.asarray(w, dtype=np.float32) for w in m.weights])

if __name__ == "__main__":
    width, seed, N, out = int(sys.argv[1]), int(sys.argv[2]), int(float(sys.argv[3])), sys.argv[4]
    depth = 16
    W = weights(width, depth, seed)
    fn = os.path.join(out, f"w{width}_s{seed}.npz")
    np.save(os.path.join(out, f"W_w{width}_s{seed}.npy"), W)
    S = np.zeros((depth, width)); Q = np.zeros((depth, width)); done = 0
    if os.path.exists(fn):
        d = np.load(fn); S, Q, done = d["S"], d["Q"], int(d["n"])
    rng = np.random.default_rng(10_000 + seed * 7919 + done)
    bs = max(4096, (1 << 22) // width)
    t0 = time.time()
    while done < N:
        x = rng.standard_normal((bs, width), dtype=np.float32)
        for l in range(depth):
            x = np.maximum(x @ W[l], 0)
            S[l] += x.sum(0, dtype=np.float64); Q[l] += np.einsum('ij,ij->j', x, x, dtype=np.float64)
        done += bs
        if done % (bs * 200) < bs or done >= N:
            np.savez(fn, S=S, Q=Q, n=done)
    print(width, seed, done, time.time() - t0)
