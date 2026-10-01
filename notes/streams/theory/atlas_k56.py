"""Central-moment atlas with the fifth- and sixth-cumulant slices the second-order closure needs (small width only).

Two passes over the SAME sample stream (same seed): pass 1 gives the exact sample means of every layer, pass 2
accumulates central moments u = z - mean (pre-activation) and v = a - mean (post-activation), so every stored
moment is the sample central moment and cumulants follow from partitions without singleton blocks.

Per layer l (z_l = a_{l-1} @ W_l, a_l = relu(z_l), a_{-1} = x ~ N(0, I)):
  pre_mean, post_mean, gate_p
  pre_cm[p]   E[u_i^p], p = 2..6 (marginal central moments)
  pre_C       E[u_i u_j]
  pre_m3      E[u_i u_j u_k]          post_m3   E[v_i v_j v_k]   (the target: all-distinct kappa3(a))
  pre_m211    E[u_i^2 u_j u_k]        (gives kappa4 (2,1,1), (2,2), (3,1), (4) slices)
  pre_m221    E[u_i^2 u_j^2 u_k]      (kappa5 (2,2,1) slice; NEW field, not in moment_atlas_np.py)
  pre_m311    E[u_i^3 u_j u_k]        (kappa5 (3,1,1) slice; NEW)
  pre_m222    E[u_i^2 u_j^2 u_k^2]    (kappa6 (2,2,2) slice; NEW)

    python atlas_k56.py --seeds 770000 --width 32 --depth 6 --n-samples 4000000 --sample-seed 11 --out DIR
"""
import argparse, os, sys, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "experiments"))
from moment_atlas_np import mlps_from_seeds


def outer_flat(x, y):
    m, n = x.shape
    return (x[:, :, None] * y[:, None, :]).reshape(m, n * n)


def run(W, N, chunk, seed):
    L, n, _ = W.shape
    # pass 1: means
    rng = np.random.default_rng(seed)
    sz = np.zeros((L, n)); sa = np.zeros((L, n)); sg = np.zeros((L, n))
    done = 0
    while done < N:
        m = min(chunk, N - done)
        a = rng.standard_normal((m, n), dtype=np.float32)
        for l in range(L):
            z = a @ W[l]; a = np.maximum(z, 0.0)
            sz[l] += z.sum(0, dtype=np.float64); sa[l] += a.sum(0, dtype=np.float64); sg[l] += (z > 0).sum(0)
        done += m
    mz, ma = sz / N, sa / N
    # pass 2: central moments on the same stream
    keys3 = ["pre_m3", "post_m3", "pre_m211", "pre_m221", "pre_m311", "pre_m222"]
    acc = {k: np.zeros((L, n, n, n)) for k in keys3}
    C = np.zeros((L, n, n)); cm = np.zeros((7, L, n))
    rng = np.random.default_rng(seed)
    done = 0
    while done < N:
        m = min(chunk, N - done)
        a = rng.standard_normal((m, n), dtype=np.float32)
        for l in range(L):
            z = a @ W[l]; a = np.maximum(z, 0.0)
            u = (z.astype(np.float64) - mz[l]); v = (a.astype(np.float64) - ma[l])
            u2 = u * u; u3 = u2 * u
            p = u2.copy()
            for q in range(2, 7):
                cm[q, l] += p.sum(0); p *= u
            C[l] += u.T @ u
            uu = outer_flat(u, u)
            acc["pre_m3"][l] += (u.T @ uu).reshape(n, n, n)
            acc["pre_m211"][l] += (u2.T @ uu).reshape(n, n, n)
            acc["pre_m311"][l] += (u3.T @ uu).reshape(n, n, n)
            u2u = outer_flat(u2, u)
            acc["pre_m221"][l] += (u2.T @ u2u).reshape(n, n, n)
            acc["pre_m222"][l] += (u2.T @ outer_flat(u2, u2)).reshape(n, n, n)
            acc["post_m3"][l] += (v.T @ outer_flat(v, v)).reshape(n, n, n)
        done += m
    out = dict(pre_mean=mz, post_mean=ma, gate_p=sg / N, pre_C=C / N, pre_cm=cm / N, n_samples=N, sample_seed=seed, weights=W)
    out.update({k: v / N for k, v in acc.items()})
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, required=True); ap.add_argument("--count", type=int, default=1)
    ap.add_argument("--width", type=int, default=32); ap.add_argument("--depth", type=int, default=6)
    ap.add_argument("--n-samples", type=int, default=4_000_000); ap.add_argument("--chunk", type=int, default=4096)
    ap.add_argument("--sample-seed", type=int, default=11); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    for i, (name, seed, W, _) in enumerate(mlps_from_seeds(a.seeds, a.count, a.width, a.depth)):
        t0 = time.time()
        res = run(W, a.n_samples, a.chunk, a.sample_seed * 1000003 + i)
        res["name"] = name
        np.savez(os.path.join(a.out, f"mlp_{i:05d}.npz"), **res)
        print(f"{name}: width {a.width} depth {a.depth} N={a.n_samples} {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
