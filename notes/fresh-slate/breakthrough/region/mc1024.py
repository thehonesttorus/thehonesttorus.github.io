"""Monte Carlo atlas at n = 1024 (bench set w1024_d16): per-layer second moments and joint third/fourth slices,
accumulated in two independent halves (A, B) so that squared norms can be noise-corrected by cross products.

Per layer l and half h, with y = z_l - mu_ref (mu_ref = W_l^T m_{l-1}, bench truth means, m_{-1} = 0) and
t = a_l - m_ref (bench truth means), accumulated as float64 sums over samples:
  s1y (n), s2y = y^T y (n,n), s1t (n), s2t = t^T t (n,n),
  s21 = (y^2)^T y (n,n)  [E y_a^2 y_b],  s31 = (y^3)^T y [E y_a^3 y_b],  s22 = (y^2)^T (y^2) [E y_a^2 y_b^2],
  s4 = sum y^4 (n), s3 = sum y^3 (n)
Saved as one npz per MLP (float32 for the (n,n) sums divided by the half count, i.e. raw moments about mu_ref).

usage: python mc1024.py MLP_INDEX N_PER_HALF [CHUNK]
"""
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "bench"))
import bench  # noqa: E402


def run(i, n_half, chunk=4096, out=None):
    S = bench.load_set("w1024_d16")
    W = bench.weights(S, i)
    L, n, _ = W.shape
    M = S["means"][i]
    mu = np.zeros((L, n))
    mu[0] = 0.0
    for l in range(1, L):
        mu[l] = M[l - 1] @ W[l].astype(np.float64)
    keys = ["s1y", "s2y", "s1t", "s2t", "s21", "s31", "s22", "s3", "s4"] + (["s21a", "s3a"] if os.environ.get("MC_A") else [])
    acc = {h: {k: [None] * L for k in keys} for h in "AB"}
    seed = 1000 + i
    t0 = time.time()
    for h, hs in zip("AB", (0, 1)):
        rng = np.random.default_rng([seed, hs])
        done = 0
        while done < n_half:
            m = min(chunk, n_half - done)
            x = rng.standard_normal((m, n), dtype=np.float32)
            for l in range(L):
                z = x @ W[l]
                y = z - mu[l].astype(np.float32)
                a = np.maximum(z, 0.0)
                t = a - M[l].astype(np.float32)
                y2 = y * y
                upd = dict(s1y=y.sum(0, dtype=np.float64), s2y=y.T @ y, s1t=t.sum(0, dtype=np.float64), s2t=t.T @ t,
                           s21=y2.T @ y, s31=(y2 * y).T @ y, s22=y2.T @ y2,
                           s3=(y2 * y).sum(0, dtype=np.float64), s4=(y2 * y2).sum(0, dtype=np.float64))
                if "s21a" in keys:
                    t2 = t * t
                    upd["s21a"] = t2.T @ t; upd["s3a"] = (t2 * t).sum(0, dtype=np.float64)
                for k, v in upd.items():
                    v = v.astype(np.float64)
                    acc[h][k][l] = v if acc[h][k][l] is None else acc[h][k][l] + v
                x = a
            done += m
            print(f"mlp {i} half {h}: {done}/{n_half}  {time.time() - t0:.0f}s", flush=True)
    res = {}
    for h in "AB":
        for k in keys:
            res[f"{h}_{k}"] = (np.stack(acc[h][k]) / n_half).astype(np.float32 if acc[h][k][0].ndim == 2 else np.float64)
    res["mu_ref"] = mu
    res["m_ref"] = M
    res["n_half"] = n_half
    out = out or os.path.join(HERE, "data", f"mc1024_mlp{i}_N{n_half}{'_a' if 'MC_A' in os.environ else ''}.npz")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    np.savez(out, **res)
    print("saved", out, f"{time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    i = int(sys.argv[1]); nh = int(sys.argv[2]); ch = int(sys.argv[3]) if len(sys.argv) > 3 else 4096
    run(i, nh, ch)
