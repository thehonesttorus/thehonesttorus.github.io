"""T2: facet births transported by face-averaged arrows (FBT), oracle test of the third cumulant.

Exact telescoping (E5 in DESIGN.md).  Fix per-neuron slopes beta_s = Cov(a_s, z_s)/Var(z_s) (the face-averaged slope of
the ReLU: on a face the arrow is linear, beta is its average) and define the birth nu_s = a~_s - beta_s * z~_s (the part of
the activation not linear in its own pre-activation: it is born at the facet z_s = 0).  Then, exactly,
    z~_l = sum_{s=0}^{l-1} nu_s P_{s->l},   nu_0 = x,   P_{s->l} = W_{s+1} D_{beta_{s+1}} W_{s+2} ... D_{beta_{l-1}} W_l.
Quadratic-birth gluing: nu_{s,k} ~= c_{s,k} (z~_{s,k}^2 - Var),  c = Cov(a_k, z~_k^2) / Var(z~_k^2)  (facet density / 2 for a
Gaussian pre-activation), with the two legs of each birth carried by the linear (Gaussian) part.  Then
    kappa3(z_{l,j}) ~= 6 sum_{s<l} sum_k (P_{s->l})_{kj} c_{s,k} K_{s,l}[k,j]^2            (one birth, two legs: tree, any depth)
with K_{s,l} = Cov(z_s, z_l) (variant 'true': measured; 'lin': P_{0->s}^T P_{0->l}).  'mem' keeps only s = l-1 (memoryless).
Compare with the Monte Carlo kappa3(z_{l,j}).
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "bench"))
from bench import weights_from_seed  # noqa: E402


def run(n, N, seed=11, L=16, chunk=50000):
    W = weights_from_seed(seed, n, L).astype(np.float64)
    W32 = W.astype(np.float32)

    def sample(rng, b):
        h = rng.standard_normal((b, n)).astype(np.float32); zs = []
        for l in range(L):
            z = h @ W32[l]; zs.append(z.astype(np.float64)); h = np.maximum(z, 0)
        return zs
    # pass 1: first and second moments, cross-layer covariances
    rng = np.random.default_rng(99 + seed); done = 0
    S1 = np.zeros((L, n)); A1 = np.zeros((L, n)); Z2 = np.zeros((L, n)); AZ = np.zeros((L, n))
    ZZ = {}
    while done < N:
        zs = sample(rng, chunk)
        for l in range(L):
            z = zs[l]; a = np.maximum(z, 0)
            S1[l] += z.sum(0); A1[l] += a.sum(0); Z2[l] += (z * z).sum(0); AZ[l] += (a * z).sum(0)
            for s in range(l + 1):
                ZZ[s, l] = ZZ.get((s, l), 0) + zs[s].T @ z
        done += chunk
    mu = S1 / N; Ea = A1 / N; var = Z2 / N - mu ** 2; caz = AZ / N - Ea * mu
    beta = caz / np.maximum(var, 1e-12)
    K = {k: v / N - np.outer(mu[k[0]], mu[k[1]]) for k, v in ZZ.items()}
    # pass 2: third cumulants of z_l and birth coefficients c = Cov(a, z~^2) / Var(z~^2)
    rng = np.random.default_rng(99 + seed); done = 0
    Z3 = np.zeros((L, n)); AZZ = np.zeros((L, n)); Z4 = np.zeros((L, n))
    while done < N:
        zs = sample(rng, chunk)
        for l in range(L):
            zc = zs[l] - mu[l]; a = np.maximum(zs[l], 0) - Ea[l]
            Z3[l] += (zc ** 3).sum(0); AZZ[l] += (a * zc * zc).sum(0); Z4[l] += (zc ** 4).sum(0)
        done += chunk
    k3 = Z3 / N; m4 = Z4 / N
    c = (AZZ / N) / np.maximum(m4 - var ** 2, 1e-12)
    # propagators
    P = {}
    for l in range(1, L):
        for s in range(l):
            if s == l - 1:
                P[s, l] = W[l]
            else:
                P[s, l] = P[s, l - 1] * beta[l - 1][None, :] @ W[l]
    # linear cross-covariances: z^lin_l = x P_{0->l} with P_{0->l} = W_0 D W_1 ... (x-space)
    Px = [W[0]]
    for l in range(1, L):
        Px.append((Px[-1] * beta[l - 1][None, :]) @ W[l])
    r = lambda v: float(np.sqrt(np.mean(v ** 2)))
    print(f"# n={n} N={N} seed={seed}")
    for l in range(1, L):
        pred = {"true": 0, "lin": 0, "mem": 0, "legs": 0}
        for s in range(0, l):          # births at every activation layer below l (x itself is the Gaussian source)
            Kt = K[s, l]
            Kl = Px[s].T @ Px[l]
            term_t = 6 * (P[s, l] * c[s][:, None] * Kt ** 2).sum(0)
            pred["true"] = pred["true"] + term_t
            pred["lin"] = pred["lin"] + 6 * (P[s, l] * c[s][:, None] * Kl ** 2).sum(0)
            if s == l - 1:
                pred["mem"] = pred["mem"] + term_t
            Kg = (K[s, s] * beta[s][None, :]) @ P[s, l]       # legs at the birth layer, transported by the face-averaged arrows
            pred["legs"] = pred["legs"] + 6 * (P[s, l] * c[s][:, None] * Kg ** 2).sum(0)
        out = f"layer {l+1:2d}: k3 rms {r(k3[l]):.3f}"
        for k, v in pred.items():
            out += f" | {k}: relerr {r(v - k3[l]) / r(k3[l]):.3f} corr {np.corrcoef(v, k3[l])[0, 1]:.3f}"
        print(out, flush=True)


if __name__ == "__main__":
    run(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 11)
