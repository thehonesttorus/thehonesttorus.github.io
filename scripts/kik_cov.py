"""Monte Carlo second-moment matrices of the post-activations at selected layers (one float32 forward pass per batch).

  python scripts/kik_cov.py DATA NET TASK NBATCH        writes $OUT/kcov_{NET}_{TASK}.npz
Accumulates, for each layer l in LAYERS: S1 = sum h_l (n,), S2 = sum h_l h_l^T (n, n) (float64 in the task, float32 on
transfer), and the pre-activation third/fourth power sums of the next layer for reference."""
import sys, os, time, numpy as np

BS, LAYERS = 4096, (0, 1, 7, 8, 13, 14)

if __name__ == "__main__":
    D, net, task, nb = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    W = [np.ascontiguousarray(w, dtype=np.float32) for w in np.load(f"{D}/W_off{net}.npy")]; n = W[0].shape[0]
    S1 = {l: np.zeros(n) for l in LAYERS}; S2 = {l: np.zeros((n, n)) for l in LAYERS}
    rng = np.random.default_rng(5_000_011 * (net + 1) + task); t0 = time.time()
    for b in range(nb):
        h = rng.standard_normal((n, BS), dtype=np.float32)
        for l, Wl in enumerate(W):
            h = np.maximum(Wl @ h, 0)
            if l in S1:
                S1[l] += h.sum(1, dtype=np.float64); S2[l] += h @ h.T
            if l >= max(LAYERS): break
    np.savez(f"{os.environ.get('OUT', '.')}/kcov_{net}_{task}.npz", N=nb * BS, layers=np.array(LAYERS),
             **{f"S1_{l}": S1[l] for l in LAYERS}, **{f"S2_{l}": S2[l].astype(np.float32) for l in LAYERS})
    print(f"net {net} task {task}: {nb * BS} samples in {time.time() - t0:.0f}s", flush=True)
