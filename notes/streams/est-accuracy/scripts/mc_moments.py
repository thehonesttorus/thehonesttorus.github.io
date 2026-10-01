"""Monte-Carlo per-neuron pre-activation cumulants (mean, var, kappa3, kappa4) at every layer of baked MLPs,
plus the post-activation means.  Two independent half-samples are kept so noise can be estimated.

  python mc_moments.py --dataset DIR --mlps 0,1 --n 400000 --out results/mc
"""
import argparse, json, os, time
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--dataset", required=True)
ap.add_argument("--mlps", default="0")
ap.add_argument("--n", type=int, default=400000)
ap.add_argument("--chunk", type=int, default=4096)
ap.add_argument("--seed", type=int, default=12345)
ap.add_argument("--out", required=True)
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
cache = os.path.join(a.dataset, "npy")
for i in [int(x) for x in a.mlps.split(",")]:
    W = np.load(os.path.join(cache, f"w{i}.npy"))  # (L, n_in, n_out), x @ w
    L, n, _ = W.shape
    rng = np.random.default_rng(a.seed + i)
    S = np.zeros((2, 5, L, n))  # half, power 0..4 raw sums of z ; plus post-act mean sum in [.,0] slot reuse
    H = np.zeros((2, L, n))
    t0 = time.time()
    done = 0
    while done < a.n:
        m = min(a.chunk, a.n - done)
        half = 0 if done < a.n // 2 else 1
        h = rng.standard_normal((m, n), dtype=np.float32)
        for l in range(L):
            z = h @ W[l]
            z64 = z.astype(np.float64)
            z2 = z64 * z64
            S[half, 1, l] += z64.sum(0)
            S[half, 2, l] += z2.sum(0)
            S[half, 3, l] += (z2 * z64).sum(0)
            S[half, 4, l] += (z2 * z2).sum(0)
            h = np.maximum(z, 0.0)
            H[half, l] += h.sum(0, dtype=np.float64)
        S[half, 0] += m
        done += m
    out = {}
    for tag, sl in (("A", [0]), ("B", [1]), ("all", [0, 1])):
        s = S[sl].sum(0)
        cnt = s[0]
        m1, m2, m3, m4 = s[1] / cnt, s[2] / cnt, s[3] / cnt, s[4] / cnt
        var = m2 - m1 ** 2
        k3 = m3 - 3 * m1 * m2 + 2 * m1 ** 3
        c4 = m4 - 4 * m1 * m3 + 6 * m1 ** 2 * m2 - 3 * m1 ** 4
        k4 = c4 - 3 * var ** 2
        out[f"mean_{tag}"] = m1; out[f"var_{tag}"] = var; out[f"k3_{tag}"] = k3; out[f"k4_{tag}"] = k4
        out[f"hmean_{tag}"] = H[sl].sum(0) / cnt
    np.savez(os.path.join(a.out, f"mc_{i}.npz"), n=a.n, **out)
    print(i, "done", time.time() - t0, flush=True)
