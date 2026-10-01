"""Local moment atlas (numpy only) for prototyping representations at small width.

For each MLP of a baked dataset directory (schema 3.0 parquet) or for seed-regenerated MLPs,
sample N standard-normal inputs in chunks and accumulate per layer (l = 0..L-1, z_l = a_{l-1} @ W_l,
a_l = relu(z_l)):
  marginals  pre_s[p,l,:] = E[z^p], post_s[p,l,:] = E[a^p], p = 1..4;  gate_p = P(z > 0)
  pairs      pre_M11 = E[z_i z_j], pre_M21 = E[z_i^2 z_j], pre_M22 = E[z_i^2 z_j^2],
             post_M11 = E[a_i a_j], post_M21 = E[a_i^2 a_j]
  gates      GG = P(z_i > 0, z_j > 0),  GX = E[1[z_i > 0] a_{l-1,j}]  (a_{-1} = x)
Conventions follow keenanpepper's community atlases so that code written here transfers to the
1024-wide atlases unchanged.  fp32 forward pass, fp64 accumulation.

    python moment_atlas_np.py --dataset ../../scratch/dev128 --split dev --n-samples 2000000 --out atlas128
    python moment_atlas_np.py --seeds 770000 --count 2 --width 128 --depth 16 --n-samples 1000000 --out atlas128s
"""
import argparse, json, os, time
import numpy as np


def mlps_from_dataset(d, split):
    import pyarrow.parquet as pq
    import glob
    files = sorted(glob.glob(os.path.join(d, "data", f"{split}-*.parquet")))
    for f in files:
        t = pq.read_table(f, columns=["mlp_name", "mlp_seed", "weights", "all_layer_means"])
        for r in range(t.num_rows):
            W = np.asarray(t.column("weights")[r].as_py(), dtype=np.float32)
            yield str(t.column("mlp_name")[r].as_py()), int(t.column("mlp_seed")[r].as_py()), W, \
                np.asarray(t.column("all_layer_means")[r].as_py(), dtype=np.float32)


def mlps_from_seeds(seed0, count, width, depth):
    for i in range(count):
        rng = np.random.default_rng(seed0 + i)
        W = (rng.standard_normal((depth, width, width)) * np.sqrt(2.0 / width)).astype(np.float32)
        yield f"seed-{seed0 + i}", seed0 + i, W, None


def atlas(W, n_samples, chunk, sample_seed, pairs=True, gates=True, k3=False, k4=False):
    L, n, _ = W.shape
    if k3:
        pre_M3 = np.zeros((L, n, n, n)); post_M3 = np.zeros((L, n, n, n))   # full raw third moments (small width only)
    if k4:
        pre_M211 = np.zeros((L, n, n, n))   # E[z_i^2 z_j z_k]: the (2,1,1) raw fourth moment of the pre-activation
    pre_s = np.zeros((4, L, n)); post_s = np.zeros((4, L, n)); gate = np.zeros((L, n))
    if pairs:
        pM11 = np.zeros((L, n, n)); pM21 = np.zeros((L, n, n)); pM22 = np.zeros((L, n, n))
        aM11 = np.zeros((L, n, n)); aM21 = np.zeros((L, n, n))
    if gates:
        GG = np.zeros((L, n, n)); GX = np.zeros((L, n, n))
    rng = np.random.default_rng(sample_seed)
    done = 0
    while done < n_samples:
        m = min(chunk, n_samples - done)
        a_prev = rng.standard_normal((m, n), dtype=np.float32)
        for l in range(L):
            z = a_prev @ W[l]
            a = np.maximum(z, 0.0)
            zd = z.astype(np.float64); ad = a.astype(np.float64)
            g = (z > 0).astype(np.float64)
            zp, ap = zd.copy(), ad.copy()
            for p in range(4):
                pre_s[p, l] += zp.sum(0); post_s[p, l] += ap.sum(0)
                zp *= zd; ap *= ad
            gate[l] += g.sum(0)
            if pairs:
                z2 = zd * zd
                pM11[l] += zd.T @ zd; pM21[l] += z2.T @ zd; pM22[l] += z2.T @ z2
                aM11[l] += ad.T @ ad; aM21[l] += (ad * ad).T @ ad
            if gates:
                GG[l] += g.T @ g
                GX[l] += g.T @ a_prev.astype(np.float64)
            if k3:
                # one GEMM per tensor: (n, m) @ (m, n*n); the (m, n, n) outer-product block is the memory limit,
                # so k3 runs use a smaller chunk (see main)
                # fp32 GEMM per chunk (chunk-level rounding ~1e-6 relative, far below the MC noise), fp64 accumulation
                pre_M3[l] += (z.T @ (z[:, :, None] * z[:, None, :]).reshape(m, n * n)).reshape(n, n, n)
                post_M3[l] += (a.T @ (a[:, :, None] * a[:, None, :]).reshape(m, n * n)).reshape(n, n, n)
            if k4:
                pre_M211[l] += ((z * z).T @ (z[:, :, None] * z[:, None, :]).reshape(m, n * n)).reshape(n, n, n)
            a_prev = a
        done += m
    N = float(n_samples)
    out = dict(pre_s=pre_s / N, post_s=post_s / N, gate_p=gate / N, n_samples=n_samples, sample_seed=sample_seed)
    if pairs:
        out.update(pre_M11=(pM11 / N).astype(np.float32), pre_M21=(pM21 / N).astype(np.float32),
                   pre_M22=(pM22 / N).astype(np.float32), post_M11=(aM11 / N).astype(np.float32),
                   post_M21=(aM21 / N).astype(np.float32))
    if gates:
        out.update(gate_GG=(GG / N).astype(np.float32), gate_GX=(GX / N).astype(np.float32))
    if k3:
        out.update(pre_M3=pre_M3 / N, post_M3=post_M3 / N)
    if k4:
        out.update(pre_M211=pre_M211 / N)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset"); ap.add_argument("--split", default="dev")
    ap.add_argument("--seeds", type=int); ap.add_argument("--count", type=int, default=1)
    ap.add_argument("--width", type=int, default=128); ap.add_argument("--depth", type=int, default=16)
    ap.add_argument("--n-samples", type=int, default=1_000_000); ap.add_argument("--chunk", type=int, default=16384)
    ap.add_argument("--no-pairs", action="store_true"); ap.add_argument("--no-gates", action="store_true")
    ap.add_argument("--k3", action="store_true", help="also accumulate the full third-moment tensors (width <= 256)")
    ap.add_argument("--k4", action="store_true", help="also accumulate the (2,1,1) fourth-moment tensor of the pre-activation")
    ap.add_argument("--limit", type=int, default=None, help="process only the first LIMIT MLPs (after --skip)")
    ap.add_argument("--skip", type=int, default=0, help="skip the first SKIP MLPs of the source")
    ap.add_argument("--sample-seed", type=int, default=20260824,
                    help="base sample seed; two runs with different seeds give independent atlases of the same MLPs")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    src = mlps_from_dataset(args.dataset, args.split) if args.dataset else mlps_from_seeds(args.seeds, args.count, args.width, args.depth)
    for i, (name, seed, W, gt) in enumerate(src):
        if i < args.skip:
            continue
        if args.limit is not None and i >= args.skip + args.limit:
            break
        t0 = time.time()
        chunk = min(args.chunk, 1024) if (args.k3 or args.k4) else args.chunk
        res = atlas(W, args.n_samples, chunk, sample_seed=args.sample_seed ^ (i * 2654435761 % 2**63),
                    pairs=not args.no_pairs, gates=not args.no_gates, k3=args.k3, k4=args.k4)
        res.update(name=name, mlp_seed=seed, weights=W)
        if gt is not None:
            res["gt_mean"] = gt
        np.savez(os.path.join(args.out, f"mlp_{i:05d}.npz"), **res)
        err = float(np.sqrt(np.mean((res["post_s"][0] - gt) ** 2))) if gt is not None else float("nan")
        print(f"{name}: N={args.n_samples} width={W.shape[1]} depth={W.shape[0]} {time.time()-t0:.0f}s  rms(mean - gt)={err:.2e}", flush=True)


if __name__ == "__main__":
    main()
