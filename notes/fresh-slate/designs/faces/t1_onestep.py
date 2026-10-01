"""T1: one-step oracle test of the barycentric face closure (BFC).

From a large Monte Carlo sample at layer l (truth), compute the exact face state (p, m, v, Gamma, X, Xi), apply one BFC
step (E4 + gluings G1, G2) and compare with the true layer-(l+1) quantities:
  (a) third / fourth cumulants of z_{l+1,j}: true vs gate-field-only (G1, true gate law) vs G2 (pairwise face closure)
  (b) one-step readout error of E a_{l+1}: Gaussian with exact variance, BFC (G2 cumulants), oracle (true cumulants).
Usage: python t1_onestep.py --n 128 --N 4000000 --layers 1,3,7,11,15 --seed 11
"""
import argparse
import os
import sys
import time

import numpy as np
from scipy.special import ndtr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "bench"))
from bench import weights_from_seed  # noqa: E402


def phi(t):
    return np.exp(-0.5 * t * t) / np.sqrt(2 * np.pi)


def relu_readout(mu, var, k3=0.0, k4=0.0):
    """E relu(z) for z with mean mu, variance var, cumulants k3, k4 (Edgeworth; the Hermite terms are facet densities)."""
    s = np.sqrt(var); t = mu / s; f = phi(t)
    base = mu * ndtr(t) + s * f
    return base - k3 * t * f / (6 * s ** 2) + k4 * (t * t - 1) * f / (24 * s ** 3) + k3 ** 2 * (t ** 4 - 6 * t * t + 3) * f / (72 * s ** 5)


def g2_cumulants(Lam, p, Gam):
    """Cumulants 3, 4 of U_j = Lam[:, j] . (g - p) from 1- and 2-body face data only (all-distinct terms dropped)."""
    q = 1 - 2 * p
    pq = p * (1 - p)
    k3g = pq * q
    k4g = pq * (1 - 6 * pq)
    G0 = Gam - np.diag(np.diag(Gam))           # off-diagonal pair covariances
    GL = G0 @ Lam                               # (n, n): sum_k G0[i,k] Lam[k,j]
    L2 = Lam ** 2
    k3 = (Lam ** 3 * k3g[:, None]).sum(0) + 3 * (L2 * q[:, None] * GL).sum(0)
    s31 = 1 - 6 * p + 6 * p * p
    M22 = np.outer(q, q) * G0 - 2 * G0 ** 2
    k4 = (Lam ** 4 * k4g[:, None]).sum(0) + 4 * (Lam ** 3 * s31[:, None] * GL).sum(0) + 3 * (L2 * (M22 @ L2)).sum(0)
    return k3, k4


class Acc:
    """Streaming accumulators of the face state at one layer and cumulants at the next."""

    def __init__(self, n):
        self.N = 0
        self.sg = np.zeros(n); self.sa = np.zeros(n); self.sa2 = np.zeros(n)
        self.gg = np.zeros((n, n)); self.ga = np.zeros((n, n)); self.aa = np.zeros((n, n))
        self.zp = [np.zeros(n) for _ in range(4)]     # raw moments of z_{l+1}
        self.ap = np.zeros(n)

    def add(self, z, z1):
        g = (z > 0).astype(np.float64); a = np.maximum(z, 0).astype(np.float64)
        self.N += len(z)
        self.sg += g.sum(0); self.sa += a.sum(0); self.sa2 += (a * a).sum(0)
        self.gg += g.T @ g; self.ga += g.T @ a; self.aa += a.T @ a
        z1 = z1.astype(np.float64)
        for k in range(4):
            self.zp[k] += (z1 ** (k + 1)).sum(0)
        self.ap += np.maximum(z1, 0).sum(0)


def run(n, N, layers, seed, chunk=50000):
    W = weights_from_seed(seed, n, 16).astype(np.float32)
    rng = np.random.default_rng(1234 + seed)
    accs = {l: Acc(n) for l in layers}
    done = 0
    t0 = time.time()
    while done < N:
        b = min(chunk, N - done)
        h = rng.standard_normal((b, n)).astype(np.float32)
        zs = []
        for l in range(16):
            z = h @ W[l]; zs.append(z); h = np.maximum(z, 0)
        for l in layers:
            accs[l].add(zs[l], zs[l + 1])
        done += b
    print(f"# n={n} N={N} seed={seed} sampled in {time.time() - t0:.0f}s", flush=True)
    rows = []
    for l in layers:
        A = accs[l]; NN = A.N
        p = A.sg / NN; Ea = A.sa / NN
        m = Ea / np.maximum(p, 1e-12)
        Gam = A.gg / NN - np.outer(p, p)
        # r_k = g_k (z_k - m_k) = a_k - m_k g_k
        Ega = A.ga / NN                                   # E[g_i a_k]
        X = Ega - np.outer(p, Ea) - (A.gg / NN - np.outer(p, p)) * m[None, :]   # Cov(g_i, a_k) - m_k Cov(g_i, g_k)
        Caa = A.aa / NN - np.outer(Ea, Ea)
        Wl = W[l + 1].astype(np.float64)
        # exact next-layer moments (from samples)
        M = [A.zp[k] / NN for k in range(4)]
        mu1 = M[0]; c2 = M[1] - mu1 ** 2
        c3 = M[2] - 3 * mu1 * M[1] + 2 * mu1 ** 3
        c4 = M[3] - 4 * mu1 * M[2] + 6 * mu1 ** 2 * M[1] - 3 * mu1 ** 4 - 3 * c2 ** 2
        Ea1 = A.ap / NN
        # BFC arrow
        mu_p = Ea @ Wl
        Cp = Wl.T @ Caa @ Wl
        reg = 1e-9 * np.eye(n)
        Lam = m[:, None] * Wl + np.linalg.solve(Gam + reg, X @ Wl)
        k3_g2, k4_g2 = g2_cumulants(Lam, p, Gam)
        var_p = np.diag(Cp)
        e_gauss = relu_readout(mu_p, var_p) - Ea1
        e_bfc3 = relu_readout(mu_p, var_p, k3_g2, 0 * k4_g2) - Ea1
        e_bfc = relu_readout(mu_p, var_p, k3_g2, k4_g2) - Ea1
        e_orc3 = relu_readout(mu1, c2, c3, 0 * c4) - Ea1
        e_orc = relu_readout(mu1, c2, c3, c4) - Ea1
        # how much of the true kappa3 is linearly explained by G2
        def rel(a, b):
            return float(np.sqrt(np.mean((a - b) ** 2) / np.mean(b ** 2)))
        row = dict(layer=l + 2, k3_rms=float(np.sqrt(np.mean(c3 ** 2))), k3_g2_relerr=rel(k3_g2, c3),
                   k3_corr=float(np.corrcoef(k3_g2, c3)[0, 1]),
                   k4_rms=float(np.sqrt(np.mean(c4 ** 2))), k4_g2_relerr=rel(k4_g2, c4),
                   mse_gauss=float(np.mean(e_gauss ** 2)), mse_bfc_k3=float(np.mean(e_bfc3 ** 2)),
                   mse_bfc=float(np.mean(e_bfc ** 2)), mse_orc_k3=float(np.mean(e_orc3 ** 2)), mse_orc=float(np.mean(e_orc ** 2)),
                   var_check=float(np.max(np.abs(var_p - c2))))
        rows.append(row)
        print(" ".join(f"{k}={v:.3e}" if isinstance(v, float) else f"{k}={v}" for k, v in row.items()), flush=True)
    return rows


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=128)
    ap.add_argument("--N", type=int, default=2000000)
    ap.add_argument("--layers", default="0,2,6,10,14")
    ap.add_argument("--seed", type=int, default=11)
    a = ap.parse_args()
    run(a.n, a.N, [int(x) for x in a.layers.split(",")], a.seed)
