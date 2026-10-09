"""Monte Carlo second-, third- and fourth-order index slices of the pre-activations z_l for selected layers:
  python scripts/mc_slices.py NET N SEED --layers 1,3,5   -> $OUT/mcs_{NET}_{SEED}.npz with, per listed layer l,
  M2 = sum z z^T, M21 = sum (z^2)^T z, M22 = sum (z^2)^T (z^2), M31 = sum (z^3)^T z, and per-neuron sums s1..s4, h1."""
import sys, os, time, numpy as np
net, N, seed = int(sys.argv[1]), int(float(sys.argv[2])), int(sys.argv[3])
layers = [int(x) for x in sys.argv[sys.argv.index("--layers") + 1].split(",")]
D = os.environ.get("DATA", "."); OUT = os.environ.get("OUT", ".")
W = [np.ascontiguousarray(w.T).astype(np.float32) for w in np.load(f"{D}/W_off{net}.npy")]
L, n = len(W), W[0].shape[0]
S = {k: np.zeros((L, n)) for k in ("s1", "s2", "s3", "s4", "h1")}
M = {l: {k: np.zeros((n, n)) for k in ("M2", "M21", "M22", "M31")} for l in layers}
rng = np.random.default_rng(seed); B = 4096; t0 = time.time()
for it in range(N // B):
    h = rng.standard_normal((B, n), dtype=np.float32)
    for l in range(L):
        z = h @ W[l]; h = np.maximum(z, 0.0)
        z64 = z.astype(np.float64); z2 = z64 * z64; S["s1"][l] += z64.sum(0); S["s2"][l] += z2.sum(0)
        S["s3"][l] += (z2 * z64).sum(0); S["s4"][l] += (z2 * z2).sum(0); S["h1"][l] += h.astype(np.float64).sum(0)
        if l in M:
            m = M[l]; m["M2"] += z64.T @ z64; m["M21"] += z2.T @ z64; m["M22"] += z2.T @ z2; m["M31"] += (z2 * z64).T @ z64
    if it % 100 == 0: print(f"net {net} seed {seed}: {(it + 1) * B} samples at {time.time() - t0:.0f}s", flush=True)
np.savez(f"{OUT}/mcs_{net}_{seed}.npz", N=(N // B) * B, layers=np.array(layers), **S,
         **{f"{k}_{l}": M[l][k] for l in layers for k in M[l]})
print(f"done {(N // B) * B} samples in {time.time() - t0:.0f}s", flush=True)
