"""ARC mean propagation (K=1) against WhestBench Phase 2 ground truth.

Reads only the small columns (seed, name, all_layer_means, avg_variance) of one
shard of the public ``mini`` split via HTTP range requests, regenerates each
network's weights from its ``mlp_seed`` with the dataset's seed protocol
(``SeedSequence(seed).spawn(3)[0]`` feeding ``default_rng``; He scale
``sqrt(2/width)``; float32 cast; see ``whestbench/seeds.py`` and
``whestbench/generation.py``), verifies the regeneration against the exact
layer-1 means ``|w_j|/sqrt(2 pi)``, and compares mean propagation with the
shipped 1e9-sample ``all_layer_means``.

Needs network access to huggingface.co and ~200 MB of memory (the 64 MB
weight tensors are regenerated, never downloaded).  Recorded output for the
first shard is in ``phase2_output.txt``.

Requires: numpy, scipy, pyarrow, fsspec, aiohttp (or requests).
"""
import time

import fsspec
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
from scipy.stats import norm

URL = ("https://huggingface.co/datasets/aicrowd/arc-whestbench-public-2026/resolve/"
       "v2-phase2/data/mini-00000-of-00007.parquet")
DEPTH, WIDTH = 16, 1024
SQ2PI = np.sqrt(2 * np.pi)
BUDGET_SAMPLES = 65536  # 2^41 FLOPs / (2 * depth * width^2) forward-pass FLOPs


def to_np(col):
    a = col.combine_chunks() if isinstance(col, pa.ChunkedArray) else col
    while pa.types.is_list(a.type) or pa.types.is_fixed_size_list(a.type) or pa.types.is_large_list(a.type):
        a = a.flatten()
    return a.to_numpy(zero_copy_only=False)


def regenerate_weights(seed):
    ss = np.random.SeedSequence(int(seed)).spawn(3)      # ss[0]: weight stream
    rng = np.random.default_rng(ss[0])
    scale = float(np.sqrt(2.0 / WIDTH))
    return [(rng.standard_normal((WIDTH, WIDTH)) * scale).astype(np.float32) for _ in range(DEPTH)]


def mean_propagation(Ws):
    """ARC Algorithm 1: track the mean vector and the mean of the covariance trace."""
    mu = np.zeros(WIDTH); s2mean = 1.0; out = []
    for W in Ws:
        Wl = W.astype(np.float64).T           # h_l = max(0, h_{l-1} @ W_l): neuron j uses column j
        mu = Wl @ mu; s2 = s2mean * np.sum(Wl ** 2, axis=1); s = np.sqrt(s2)
        a = mu * norm.cdf(mu / s) + s * norm.pdf(mu / s)                   # E[phi(Y)], Y ~ N(mu, s^2)
        b = (mu ** 2 + s2) * norm.cdf(mu / s) + mu * s * norm.pdf(mu / s)  # E[phi(Y)^2]
        out.append(a.copy()); mu = a; s2mean = np.mean(b - a ** 2)
    return np.array(out)


def main():
    t0 = time.time()
    fs = fsspec.filesystem("https")
    with fs.open(URL, block_size=2 ** 20) as f:
        pf = pq.ParquetFile(f)
        tb = pf.read_row_group(0, columns=["mlp_name", "mlp_seed", "all_layer_means", "avg_variance"])
    Y = to_np(tb.column("all_layer_means")).astype(np.float64).reshape(tb.num_rows, DEPTH, WIDTH)
    names = tb.column("mlp_name").to_pylist(); seeds = tb.column("mlp_seed").to_pylist()
    avs = tb.column("avg_variance").to_pylist()
    print(f"read {tb.num_rows} rows in {time.time() - t0:.1f}s")
    M, Ym = [], []
    for k in range(tb.num_rows):
        Ws = regenerate_weights(seeds[k])
        l1 = np.linalg.norm(Ws[0].astype(np.float64), axis=0) / SQ2PI   # exact layer-1 means
        P = mean_propagation(Ws); mse = np.mean((P - Y[k]) ** 2, axis=1)
        M.append(mse); Ym.append(np.mean(Y[k] ** 2, axis=1))
        print(f"mlp {k:2d} {names[k]:22s} seed={seeds[k]}  layer-1 regen check max|diff|={np.abs(l1 - Y[k, 0]).max():.1e}"
              f" | final MSE={mse[-1]:.3e} var={avs[k]:.3e} MSE/var={mse[-1] / avs[k]:.2e}"
              f" = {mse[-1] / (avs[k] / BUDGET_SAMPLES):.0f}x budget-matched MC")
    M = np.array(M); Ym = np.array(Ym)
    print(f"\nper-layer, averaged over {len(M)} networks: mean propagation K=1")
    for l in range(DEPTH):
        print(f"  layer {l + 1:2d}: MSE={M[:, l].mean():.3e}  MSE/mean^2={np.mean(M[:, l] / Ym[:, l]):.2e}")
    print(f"total {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
