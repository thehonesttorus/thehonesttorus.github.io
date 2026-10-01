"""T2: Gaussian closure variants at n=1024: linearised cross-covariance (K=1) vs Gaussian-exact bivariate ReLU covariance
truncated at Hermite order K (Cov_pq = sum_k s_p s_q h_k(a_p) h_k(a_q) rho^k / k!), with the exact chi radial factor."""
import sys, json, time
import numpy as np
from common import bench, phi, Phi, chi_mean_ratio

def hk(a, K):
    p = phi(a); H = [Phi(a), p]
    he = [np.ones_like(a), a]          # probabilists' Hermite He_0, He_1
    for k in range(3, K + 1):
        j = k - 2                       # h_k = (-1)^j He_j(a) phi(a)
        while len(he) <= j: he.append(a * he[-1] - (len(he) - 1) * he[-2])
        H.append(((-1) ** j) * he[j] * p)
    return H[:K]

def closure(W, K=1, chi=False):
    W = W.astype(np.float64); L, n, _ = W.shape
    out = []; m = np.zeros(n); S = W[0].T @ W[0]
    fact = [1.0]
    for k in range(1, K + 1): fact.append(fact[-1] * k)
    for l in range(L):
        if l > 0: m = mu @ W[l]; S = W[l].T @ C @ W[l]
        v = np.maximum(np.diag(S), 1e-300); s = np.sqrt(v); a = m / s
        P = Phi(a); p = phi(a)
        mu = m * P + s * p; sec = (m * m + v) * P + m * s * p
        out.append(mu)
        R = S / s[:, None] / s[None, :]
        H = hk(a, K)
        C = np.zeros_like(S); Rk = np.ones_like(S)
        for k in range(1, K + 1):
            Rk = Rk * R
            C += np.outer(s * H[k - 1], s * H[k - 1]) * Rk / fact[k]
        np.fill_diagonal(C, np.maximum(sec - mu * mu, 0.0))
    out = np.stack(out)
    return out * (chi_mean_ratio(n) if chi else 1.0)

if __name__ == "__main__":
    name = sys.argv[1]; Ks = [int(k) for k in sys.argv[2].split(",")]
    S = bench.load_set(name); res = {}
    for i in range(len(S["seeds"])):
        W = bench.weights(S, i); T = S["means"][i]; nz = S["noise"][i]
        row = []
        for K in Ks:
            for chi in (False, True):
                g = closure(W, K, chi)
                e = g[-1] - T[-1]
                res.setdefault((K, chi), []).append(((e ** 2).mean() - nz, e.mean()))
                row.append(f"K{K}{'c' if chi else ''} {((e**2).mean()-nz):.3e} (bias {e.mean():+.1e})")
        print(i, " | ".join(row), flush=True)
    for k, v in res.items():
        v = np.array(v); print(k, f"raw {v[:,0].mean():.3e}  bias {v[:,1].mean():+.2e}")
