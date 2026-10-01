"""T3: exact age attribution of the skew.  With beta_s = Cov(a_s,z_s)/Var(z_s) and births nu_s = a~_s - beta_s z~_s,
z~_l = x P_{x->l} + sum_{s<l} nu_s P_{s->l} exactly, so kappa3(z_{l,j}) = sum_{s} T_s[j],  T_s = E[(nu_s P_{s->l})_j z~_{l,j}^2]
(+ the x term).  T_s is the part of the final-layer third cumulant born at the facets of layer s and still present at l
(net of everything that later layers did to it).  Also reported: the 'self' part kappa3(nu_s P) (what a tree transport
would carry) to show how much later facets fold away."""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "bench"))
from bench import weights_from_seed

def run(n, N, seed, tgt, L=16, chunk=50000):
    W = weights_from_seed(seed, n, L).astype(np.float64); W32 = W.astype(np.float32)
    def sample(rng, b):
        x = rng.standard_normal((b, n)).astype(np.float32); h = x; zs = []
        for l in range(L):
            z = h @ W32[l]; zs.append(z.astype(np.float64)); h = np.maximum(z, 0)
        return x.astype(np.float64), zs
    rng = np.random.default_rng(3 + seed); done = 0
    S1 = np.zeros((L, n)); A1 = np.zeros((L, n)); Z2 = np.zeros((L, n)); AZ = np.zeros((L, n))
    while done < N:
        x, zs = sample(rng, chunk)
        for l in range(L):
            z = zs[l]; a = np.maximum(z, 0)
            S1[l] += z.sum(0); A1[l] += a.sum(0); Z2[l] += (z * z).sum(0); AZ[l] += (a * z).sum(0)
        done += chunk
    mu = S1 / N; Ea = A1 / N; var = Z2 / N - mu ** 2; beta = (AZ / N - Ea * mu) / np.maximum(var, 1e-12)
    for l in tgt:
        P = {}
        P[l - 1] = W[l]
        for s in range(l - 2, -1, -1):
            P[s] = (W[s + 1] * beta[s][:, None].T) if False else None
        # P_{s->l} = W_{s+1} D_{beta_{s+1}} ... W_l  (index convention: nu_s lives at activation index s, z index l)
        P = {l - 1: W[l]}
        for s in range(l - 2, -1, -1):
            P[s] = W[s + 1] @ (beta[s + 1][:, None] * P[s + 1])
        Px = W[0] @ (beta[0][:, None] * P[0])
        rng = np.random.default_rng(3 + seed); done = 0
        T = np.zeros((l + 1, n)); Sf = np.zeros((l + 1, n)); k3 = np.zeros(n); chk = 0
        while done < N:
            x, zs = sample(rng, chunk)
            zl = zs[l] - mu[l]; zl2 = zl * zl
            ys = [x @ Px]
            for s in range(l):
                nu = np.maximum(zs[s], 0) - Ea[s] - beta[s] * (zs[s] - mu[s])
                ys.append(nu @ P[s])
            for t, y in enumerate(ys):
                T[t] += (y * zl2).sum(0); Sf[t] += (y ** 3).sum(0)
            k3 += (zl ** 3).sum(0)
            chk = max(chk, np.abs(sum(ys) - zl).max())
            done += chunk
        T /= N; Sf /= N; k3 /= N
        r = lambda v: float(np.sqrt(np.mean(v ** 2)))
        print(f"# n={n} seed={seed} target z-layer {l+1}: kappa3 rms {r(k3):.3f}; telescoping max err {chk:.1e}; sum_s T_s - k3 rms {r(T.sum(0)-k3):.1e}")
        print("  age(l-s) | T_s rms | T_s . k3 / |k3|^2 (share) | self kappa3(nu_s P) rms")
        for t in range(l + 1):
            age = l - (t - 1) if t > 0 else l + 1
            print(f"  {'x' if t == 0 else age:>8} | {r(T[t]):.4f} | {float(T[t] @ k3 / (k3 @ k3)):+.3f} | {r(Sf[t]):.4f}")

if __name__ == "__main__":
    run(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), [int(v) for v in sys.argv[4].split(",")])
