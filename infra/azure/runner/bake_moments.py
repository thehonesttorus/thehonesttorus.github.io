"""Moment atlas for one or more MLPs on a GPU: per-layer pre- and post-activation marginal raw
moments E[z^p], E[a^p] (p = 1..6), gate probabilities P(z>0), and (optionally) dense pairwise
blocks E[z_i z_j], E[z_i^2 z_j], E[z_i^2 z_j^2], E[a_i a_j], E[a_i^2 a_j] per layer; fp32 forward
pass (TF32 off, like the grader), fp64 accumulation.  Same conventions as the community atlases
(keenanpepper/*): z_l = a_{l-1} @ W_l, a_l = relu(z_l), a_{-1} = x ~ N(0, I), l = 0..15.

    python bake_moments.py --weights-from seeds --seed0 770000 --idx 0 --n-samples 100000000 \
        --pairs --out <dir or blob container URL> [--chunk 131072] [--marginal-order 6]
    python bake_moments.py --weights-from parquet --parquet shard.parquet --row 3 ...

Weights from seeds reproduce keenanpepper's d8b corpus: W = torch.randn(16,1024,1024,
generator=cpu_gen(seed0+idx)) * sqrt(2/1024).  Output: <out>/<tag>/marg.npz (+ pairs_LL.npz).
"""
import argparse, io, os, sys, time
import numpy as np
import torch


def weights_from_seed(seed, depth=16, width=1024):
    g = torch.Generator(device="cpu").manual_seed(seed)
    W = torch.randn(depth, width, width, generator=g, dtype=torch.float32)
    return W * (2.0 / width) ** 0.5


def weights_from_parquet(path, row):
    import pyarrow.parquet as pq
    t = pq.ParquetFile(path).read_row_group(0) if False else pq.read_table(path, columns=["weights", "mlp_seed", "mlp_name"])
    w = np.asarray(t.column("weights")[row].as_py(), dtype=np.float32)  # [depth, width, width] (input, output)
    return torch.from_numpy(w), int(t.column("mlp_seed")[row].as_py()), str(t.column("mlp_name")[row].as_py())


def put_bytes(dest, name, data):
    if dest.startswith("http"):
        from azure.storage.blob import ContainerClient
        ContainerClient.from_container_url(dest).upload_blob(name, data, overwrite=True)
    else:
        p = os.path.join(dest, name)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        with open(p, "wb") as f:
            f.write(data)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--weights-from", choices=["seeds", "parquet"], default="seeds")
    ap.add_argument("--seed0", type=int, default=770000)
    ap.add_argument("--idx", type=int, default=0)
    ap.add_argument("--parquet"); ap.add_argument("--row", type=int, default=0)
    ap.add_argument("--n-samples", type=int, default=100_000_000)
    ap.add_argument("--chunk", type=int, default=131072)
    ap.add_argument("--marginal-order", type=int, default=6)
    ap.add_argument("--pairs", action="store_true", help="also accumulate dense n x n pair blocks (5 per layer, 20 MB each)")
    ap.add_argument("--sample-seed", type=int, default=None)
    ap.add_argument("--out", required=True)
    ap.add_argument("--tag", default=None)
    args = ap.parse_args()
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    if args.weights_from == "seeds":
        W = weights_from_seed(args.seed0 + args.idx); mlp_seed = args.seed0 + args.idx; name = f"d8b-{args.idx:05d}"
    else:
        W, mlp_seed, name = weights_from_parquet(args.parquet, args.row)
    tag = args.tag or name
    W = W.to(dev)
    depth, width = W.shape[0], W.shape[1]
    P = args.marginal_order
    pre_s = torch.zeros(P, depth, width, dtype=torch.float64, device=dev)
    post_s = torch.zeros(P, depth, width, dtype=torch.float64, device=dev)
    gate = torch.zeros(depth, width, dtype=torch.float64, device=dev)
    if args.pairs:
        pre_M11 = torch.zeros(depth, width, width, dtype=torch.float64, device=dev)
        pre_M21 = torch.zeros_like(pre_M11); pre_M22 = torch.zeros_like(pre_M11)
        post_M11 = torch.zeros_like(pre_M11); post_M21 = torch.zeros_like(pre_M11)
    sample_seed = args.sample_seed if args.sample_seed is not None else (20260824 ^ ((args.idx * 2654435761) % 2**63))
    g = torch.Generator(device=dev).manual_seed(sample_seed)
    n_done, t0 = 0, time.time()
    while n_done < args.n_samples:
        m = min(args.chunk, args.n_samples - n_done)
        a = torch.randn(m, width, generator=g, device=dev, dtype=torch.float32)
        for l in range(depth):
            z = a @ W[l]
            a = torch.relu(z)
            zd = z.double(); ad = a.double()
            zp = zd.clone(); apw = ad.clone()
            for p in range(P):
                pre_s[p, l] += zp.sum(0); post_s[p, l] += apw.sum(0)
                zp = zp * zd; apw = apw * ad
            gate[l] += (z > 0).double().sum(0)
            if args.pairs:
                pre_M11[l] += zd.T @ zd
                z2 = zd * zd
                pre_M21[l] += z2.T @ zd
                pre_M22[l] += z2.T @ z2
                post_M11[l] += ad.T @ ad
                post_M21[l] += (ad * ad).T @ ad
        n_done += m
        if (n_done // args.chunk) % 50 == 0:
            print(f"{tag}: {n_done}/{args.n_samples} ({time.time()-t0:.0f}s)", flush=True)
    N = float(args.n_samples)
    marg = dict(pre_s=(pre_s / N).cpu().numpy(), post_s=(post_s / N).cpu().numpy(), gate_p=(gate / N).cpu().numpy(),
                n_samples=args.n_samples, sample_seed=sample_seed, mlp_seed=mlp_seed, name=name, chunk=args.chunk,
                width=width, depth=depth, torch_version=torch.__version__, seconds=time.time() - t0)
    buf = io.BytesIO(); np.savez(buf, **marg); put_bytes(args.out, f"{tag}/marg.npz", buf.getvalue())
    if args.pairs:
        for l in range(depth):
            buf = io.BytesIO()
            np.savez(buf, pre_M11=(pre_M11[l] / N).float().cpu().numpy(), pre_M21=(pre_M21[l] / N).float().cpu().numpy(),
                     pre_M22=(pre_M22[l] / N).float().cpu().numpy(), post_M11=(post_M11[l] / N).float().cpu().numpy(),
                     post_M21=(post_M21[l] / N).float().cpu().numpy())
            put_bytes(args.out, f"{tag}/pairs_{l:02d}.npz", buf.getvalue())
    print(f"done {tag} in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
