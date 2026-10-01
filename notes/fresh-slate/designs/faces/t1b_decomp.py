"""T1b: where does the third cumulant of z_{l+1,j} live?  Decompose kappa3(z'_j) = sum_{ikm} W_ij W_kj W_mj kappa(a_i,a_k,a_m)
into 1-index, 2-index and all-distinct parts (activation representation), and the same for the gate field U = Lam^T (g-p)
(face representation), with the face/residual split a = m g + r.  Truth by Monte Carlo (two-pass: means first)."""
import sys, os, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "bench"))
from bench import weights_from_seed

def run(n, N, layers, seed, chunk=50000):
    W = weights_from_seed(seed, n, 16)
    def sample(rng, b):
        h = rng.standard_normal((b, n)).astype(np.float32); zs = []
        for l in range(16):
            z = h @ W[l]; zs.append(z); h = np.maximum(z, 0)
        return zs
    # pass 1: means
    rng = np.random.default_rng(7 + seed); sa = {l: 0 for l in layers}; sg = {l: 0 for l in layers}; done = 0
    sga = {l: 0 for l in layers}; sgg = {l: 0 for l in layers}
    while done < N:
        zs = sample(rng, chunk)
        for l in layers:
            g = (zs[l] > 0).astype(np.float64); a = np.maximum(zs[l], 0).astype(np.float64)
            sa[l] = sa[l] + a.sum(0); sg[l] = sg[l] + g.sum(0); sga[l] = sga[l] + g.T @ a; sgg[l] = sgg[l] + g.T @ g
        done += chunk
    st = {}
    for l in layers:
        p = sg[l] / N; Ea = sa[l] / N; m = Ea / p
        Gam = sgg[l] / N - np.outer(p, p)
        X = sga[l] / N - np.outer(p, Ea) - Gam * m[None, :]
        Wl = W[l + 1].astype(np.float64)
        Lam = m[:, None] * Wl + np.linalg.solve(Gam + 1e-9 * np.eye(n), X @ Wl)
        st[l] = (p, Ea, m, Lam, Wl)
    # pass 2: centred third moments
    rng = np.random.default_rng(7 + seed); done = 0
    acc = {l: dict(a21=0, a3=0, z3=0, u3=0, u21=0, g3=0) for l in layers}
    while done < N:
        zs = sample(rng, chunk)
        for l in layers:
            p, Ea, m, Lam, Wl = st[l]
            a = np.maximum(zs[l], 0).astype(np.float64) - Ea; g = (zs[l] > 0) - p
            A = acc[l]
            A["a21"] = A["a21"] + (a * a).T @ a          # E[a_i^2 a_k]
            A["a3"] = A["a3"] + (a ** 3).sum(0)
            A["z3"] = A["z3"] + ((a @ Wl) ** 3).sum(0)
            A["u3"] = A["u3"] + ((g @ Lam) ** 3).sum(0)
            A["u21"] = A["u21"] + (g * g).T @ g
        done += chunk
    for l in layers:
        p, Ea, m, Lam, Wl = st[l]; A = {k: v / N for k, v in acc[l].items()}
        tot = A["z3"]
        one = ((Wl ** 3) * A["a3"][:, None]).sum(0)
        M = A["a21"].copy(); np.fill_diagonal(M, 0)
        two = 3 * ((Wl ** 2) * (M @ Wl)).sum(0)
        dist = tot - one - two
        utot = A["u3"]
        q = 1 - 2 * p; pq = p * (1 - p)
        uone = ((Lam ** 3) * (pq * q)[:, None]).sum(0)
        M = A["u21"].copy(); np.fill_diagonal(M, 0)
        utwo = 3 * ((Lam ** 2) * (M @ Lam)).sum(0)
        udist = utot - uone - utwo
        r = lambda x: float(np.sqrt(np.mean(x ** 2)))
        print(f"layer {l+2}: k3(z') rms {r(tot):.3f} | act-rep: 1-idx {r(one):.3f} 2-idx {r(two):.3f} distinct {r(dist):.3f}"
              f" | gate field U: k3 rms {r(utot):.3f} 1-idx {r(uone):.3f} 2-idx {r(utwo):.3f} distinct {r(udist):.3f}"
              f" corr(k3U,k3z) {np.corrcoef(utot, tot)[0,1]:.2f}", flush=True)

if __name__ == "__main__":
    n = int(sys.argv[1]); N = int(sys.argv[2]); layers = [int(x) for x in sys.argv[3].split(",")]
    run(n, N, layers, 11)
