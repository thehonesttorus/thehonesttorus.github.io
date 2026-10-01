"""T5: the folding of old facet births, measured.  Exact: kappa3(z_{l,j}) = T_x[j] + sum_s sum_k P_{s->l}[k,j] H_s[k,j],
H_s[k,j] = E[nu_{s,k} z~_{l,j}^2]  (the response of the final variance to the birth at facet (s,k)).
Tree value with renormalised legs: H^tree_s[k,j] = 2 c_{s,k} K_s[k,j]^2,  K_s = C_s D_{beta_s} P_{s->l},  c = E[nu z~^2]/(2 Var^2).
Per age: alpha = best scalar H ~ alpha H^tree (on the P-weighted contributions T_s[j]), R^2, and the share of T_s."""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "bench"))
from bench import weights_from_seed

def run(n, N, seed, l, L=16, chunk=25000):
    W = weights_from_seed(seed, n, L).astype(np.float64); W32 = W.astype(np.float32)
    def sample(rng, b):
        x = rng.standard_normal((b, n)).astype(np.float32); h = x; zs = []
        for i in range(L):
            z = h @ W32[i]; zs.append(z.astype(np.float64)); h = np.maximum(z, 0)
        return x.astype(np.float64), zs
    rng = np.random.default_rng(21 + seed); done = 0
    S1 = np.zeros((L, n)); A1 = np.zeros((L, n)); AZ = np.zeros((L, n)); ZZ = np.zeros((L, n, n))
    while done < N:
        x, zs = sample(rng, chunk)
        for i in range(l + 1):
            z = zs[i]; a = np.maximum(z, 0)
            S1[i] += z.sum(0); A1[i] += a.sum(0); AZ[i] += (a * z).sum(0); ZZ[i] += z.T @ z
        done += chunk
    mu = S1 / N; Ea = A1 / N
    C = ZZ / N - mu[:, :, None] * mu[:, None, :]
    var = np.einsum("lii->li", C)
    beta = (AZ / N - Ea * mu) / var
    P = {l - 1: W[l]}
    for s in range(l - 2, -1, -1):
        P[s] = W[s + 1] @ (beta[s + 1][:, None] * P[s + 1])
    rng = np.random.default_rng(21 + seed); done = 0
    H = np.zeros((l, n, n)); NZ = np.zeros((l, n)); k3 = np.zeros(n)
    while done < N:
        x, zs = sample(rng, chunk)
        zl2 = (zs[l] - mu[l]) ** 2
        for s in range(l):
            zc = zs[s] - mu[s]
            nu = np.maximum(zs[s], 0) - Ea[s] - beta[s] * zc
            H[s] += nu.T @ zl2; NZ[s] += (nu * zc * zc).sum(0)
        k3 += ((zs[l] - mu[l]) ** 3).sum(0)
        done += chunk
    H /= N; NZ /= N; k3 /= N
    r = lambda v: float(np.sqrt(np.mean(v ** 2)))
    print(f"# n={n} seed={seed} z-layer {l+1}: kappa3 rms {r(k3):.4f}")
    print(" age | share of k3 | alpha (renorm. legs) | R^2 | alpha (one-step legs, s=l-1 only)")
    for s in range(l):
        c = NZ[s] / (2 * var[s] ** 2)
        K = (C[s] * beta[s][None, :]) @ P[s]
        Ht = 2 * c[:, None] * K ** 2
        T = (P[s] * H[s]).sum(0); Tt = (P[s] * Ht).sum(0)
        alpha = float(T @ Tt / (Tt @ Tt)); R2 = 1 - float(np.sum((T - alpha * Tt) ** 2) / np.sum(T ** 2))
        print(f" {l - s:3d} | {float(T @ k3 / (k3 @ k3)):+.3f} | {alpha:+.3f} | {R2:.3f}")

if __name__ == "__main__":
    run(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
