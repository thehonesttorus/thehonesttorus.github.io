"""T1c: hub-face gluing of the all-distinct third cumulant.

kappa3(z'_j) = 1-index + 2-index + all-distinct (activation representation, as in T1b).  The all-distinct part needs
3-body data.  Gluing rule (conditional expectation onto the hub's local sigma-algebra): for a triple (i,k,m) with hub k,
E[a~_i a~_k a~_m] ~= E[ a~_k  E[a~_i | F_k] E[a~_m | F_k] ], summed over the three hub choices, with
  - 'lin'  : F_k = sigma(z_k), E[.|F_k] linear in z_k                    (barycentre of the whole line)
  - 'face' : F_k = sigma(z_k), E[.|F_k] linear on each half-face of k     (span{z_k, a_k}: the face-wise barycentric projection)
All regression coefficients use only pair data (Cov(a,z), Cov(a,a)); hub moments are 1-body.
Truth by Monte Carlo (two passes)."""
import sys, os
import numpy as np
SPHERE = "--sphere" in sys.argv
if SPHERE: sys.argv.remove("--sphere")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "bench"))
from bench import weights_from_seed

def run(n, N, layers, seed, chunk=50000):
    W = weights_from_seed(seed, n, 16).astype(np.float64)
    W32 = W.astype(np.float32)
    def sample(rng, b):
        h = rng.standard_normal((b, n)).astype(np.float32); zs = []
        if SPHERE:
            h *= np.sqrt(n) / np.linalg.norm(h, axis=1, keepdims=True)
        for l in range(16):
            z = h @ W32[l]; zs.append(z); h = np.maximum(z, 0)
        return zs
    S = {l: dict(a=0, z=0, aa=0, az=0, zz=0) for l in layers}
    rng = np.random.default_rng(5 + seed); done = 0
    while done < N:
        zs = sample(rng, chunk)
        for l in layers:
            z = zs[l].astype(np.float64); a = np.maximum(z, 0); s = S[l]
            s["a"] += a.sum(0); s["z"] += z.sum(0); s["aa"] += a.T @ a; s["az"] += a.T @ z; s["zz"] += z.T @ z
        done += chunk
    P = {}
    for l in layers:
        s = {k: v / N for k, v in S[l].items()}
        Ea, Ez = s["a"], s["z"]
        Caa = s["aa"] - np.outer(Ea, Ea); Caz = s["az"] - np.outer(Ea, Ez); Czz = s["zz"] - np.outer(Ez, Ez)
        P[l] = (Ea, Ez, Caa, Caz, Czz)
    acc = {l: dict(a21=0, a3=0, z3=0, h=0) for l in layers}
    rng = np.random.default_rng(5 + seed); done = 0
    while done < N:
        zs = sample(rng, chunk)
        for l in layers:
            Ea, Ez, Caa, Caz, Czz = P[l]
            z = zs[l].astype(np.float64); a = np.maximum(z, 0) - Ea; zc = z - Ez
            A = acc[l]
            A["a21"] = A["a21"] + (a * a).T @ a
            A["a3"] = A["a3"] + (a ** 3).sum(0)
            A["z3"] = A["z3"] + ((a @ W[l + 1]) ** 3).sum(0)
            # hub moments: E[a~_k zc_k^2], E[a~_k zc_k a~_k], E[a~_k^3]
            A["h"] = A["h"] + np.stack([(a * zc * zc).sum(0), (a * zc * a).sum(0), (a ** 3).sum(0)])
        done += chunk
    for l in layers:
        Ea, Ez, Caa, Caz, Czz = P[l]; A = {k: v / N for k, v in acc[l].items()}; Wl = W[l + 1]
        tot = A["z3"]
        one = ((Wl ** 3) * A["a3"][:, None]).sum(0)
        M = A["a21"].copy(); np.fill_diagonal(M, 0)
        two = 3 * ((Wl ** 2) * (M @ Wl)).sum(0)
        dist = tot - one - two
        vz = np.maximum(np.diag(Czz), 1e-12); va = np.diag(Caa); cza = np.diag(Caz)   # Cov(a_k, z_k)
        c_zz, c_za, c_aa = A["h"]
        preds = {}
        # lin: E[a~_i | z_k] = R_ik zc_k, R_ik = Cov(a_i, z_k)/Var z_k
        R = Caz / vz[None, :]
        np.fill_diagonal(R, 0)
        U = Wl.T @ R            # U[j,k] = sum_i W_ij R_ik  (i != k)
        D = (Wl ** 2).T @ (R ** 2)   # sum_i W_ij^2 R_ik^2
        preds["lin"] = 3 * (Wl.T * c_zz[None, :] * (U ** 2 - D)).sum(1)
        # face: E[a~_i | F_k] = Ra_ik zc_k + Rb_ik a~_k  (2x2 regression per hub k)
        G11, G12, G22 = vz, cza, va
        det = np.maximum(G11 * G22 - G12 ** 2, 1e-12)
        Ci1 = Caz; Ci2 = Caa        # Cov(a_i, z_k), Cov(a_i, a_k)
        Ra = (Ci1 * G22[None, :] - Ci2 * G12[None, :]) / det[None, :]
        Rb = (Ci2 * G11[None, :] - Ci1 * G12[None, :]) / det[None, :]
        np.fill_diagonal(Ra, 0); np.fill_diagonal(Rb, 0)
        Ua = Wl.T @ Ra; Ub = Wl.T @ Rb
        Daa = (Wl ** 2).T @ (Ra * Ra); Dab = (Wl ** 2).T @ (Ra * Rb); Dbb = (Wl ** 2).T @ (Rb * Rb)
        WT = Wl.T
        preds["face"] = 3 * (WT * (c_zz[None, :] * (Ua ** 2 - Daa) + 2 * c_za[None, :] * (Ua * Ub - Dab)
                                    + c_aa[None, :] * (Ub ** 2 - Dbb))).sum(1)
        r = lambda x: float(np.sqrt(np.mean(x ** 2)))
        out = f"layer {l+2}: k3 rms {r(tot):.3f} 1idx {r(one):.3f} 2idx {r(two):.3f} distinct {r(dist):.3f}"
        for k, v in preds.items():
            out += f" | {k}: rms {r(v):.3f} corr {np.corrcoef(v, dist)[0,1]:.2f} relerr {r(v - dist)/r(dist):.2f}"
        print(out, flush=True)

if __name__ == "__main__":
    n = int(sys.argv[1]); N = int(sys.argv[2]); layers = [int(x) for x in sys.argv[3].split(",")]
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 11
    run(n, N, layers, seed)
