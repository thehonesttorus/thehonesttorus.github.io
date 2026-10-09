"""High-precision Monte Carlo of per-neuron pre-activation moments for one official network.
  python scripts/mc_cumulants.py NET N SEED [--cov L1,L2]   -> $OUT/mcc_{NET}_{SEED}.npz with sums s1..s4 (L, n), gate sums,
  post-activation sums h1, h2, and optionally E[z z^T] sums for the listed layers. float32 forward passes, float64 sums."""
import sys, os, time, numpy as np
net, N, seed = int(sys.argv[1]), int(float(sys.argv[2])), int(sys.argv[3])
covL = [int(x) for x in sys.argv[sys.argv.index("--cov") + 1].split(",")] if "--cov" in sys.argv else []
D = os.environ.get("DATA", "."); OUT = os.environ.get("OUT", ".")
W = [np.ascontiguousarray(w.T).astype(np.float32) for w in np.load(f"{D}/W_off{net}.npy")]   # row convention: z = h @ W^T
L, n = len(W), W[0].shape[0]
S = {k: np.zeros((L, n)) for k in ("s1", "s2", "s3", "s4", "pos", "h1", "h2")}
Cz = {l: np.zeros((n, n)) for l in covL}
rng = np.random.default_rng(seed); B = 8192; t0 = time.time()
for it in range(N // B):
    h = rng.standard_normal((B, n), dtype=np.float32)
    for l in range(L):
        z = h @ W[l]; h = np.maximum(z, 0.0)
        z64 = z.astype(np.float64); S["s1"][l] += z64.sum(0); z2 = z64 * z64; S["s2"][l] += z2.sum(0)
        S["s3"][l] += (z2 * z64).sum(0); S["s4"][l] += (z2 * z2).sum(0); S["pos"][l] += (z > 0).sum(0)
        h64 = h.astype(np.float64); S["h1"][l] += h64.sum(0); S["h2"][l] += (h64 * h64).sum(0)
        if l in Cz: Cz[l] += z.T.astype(np.float64) @ z.astype(np.float64)
    if it % 50 == 0: print(f"net {net} seed {seed}: {(it + 1) * B} samples at {time.time() - t0:.0f}s", flush=True)
np.savez(f"{OUT}/mcc_{net}_{seed}.npz", N=(N // B) * B, **S, **{f"cov{l}": Cz[l] for l in covL})
print(f"done {(N // B) * B} samples in {time.time() - t0:.0f}s", flush=True)
