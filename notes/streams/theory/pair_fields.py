"""Same-sample pair fields for the leaf-resummed closure: Gamma_ij = Cov(1[z_i > 0], a_j) and the density of z_i at 0.

Regenerates the exact sample stream of an existing atlas (moment_atlas_np.py: chunk 1024 for --k3/--k4 runs, seed
stored in the npz; atlas_k56.py: chunk 4096) and accumulates, per layer,
    EGa[i, j] = E[1[z_i > 0] a_j],  Eg = P(z > 0) (checked against the atlas gate_p: same samples),  Ea = E[a],
    dens[h][i] = E[phi((z_i)/(h sigma_i)) / (h sigma_i)]   (Gaussian-kernel density of z_i at 0, bandwidths h)
Output: <out>.npz with Gam = EGa - Eg Ea^T, w2hat = dens (the dressed vertex weight E[relu''(z)] = p_z(0)).

    python pair_fields.py ATLAS.npz OUT.npz [--k56]
"""
import sys
import numpy as np


def main(path, out, k56=False):
    z = np.load(path)
    W = z["weights"].astype(np.float32)
    L, n, _ = W.shape
    N = int(z["n_samples"])
    seed = int(z["sample_seed"])
    if k56:
        chunk = 4096
        sig = np.sqrt(np.stack([np.diag(z["pre_C"][l]) for l in range(L)]))
        gate_ref = z["gate_p"]
    else:
        chunk = 1024
        sig = np.sqrt(z["pre_s"][1] - z["pre_s"][0] ** 2)
        gate_ref = z["gate_p"]
    hs = (0.02, 0.05, 0.1)
    EGa = np.zeros((L, n, n)); Eg = np.zeros((L, n)); Ea = np.zeros((L, n)); dens = np.zeros((len(hs), L, n))
    rng = np.random.default_rng(seed)
    done = 0
    while done < N:
        m = min(chunk, N - done)
        a = rng.standard_normal((m, n), dtype=np.float32)
        for l in range(L):
            zz = a @ W[l]
            a = np.maximum(zz, 0.0)
            g = (zz > 0).astype(np.float64); ad = a.astype(np.float64); zd = zz.astype(np.float64)
            EGa[l] += g.T @ ad; Eg[l] += g.sum(0); Ea[l] += ad.sum(0)
            for t, h in enumerate(hs):
                bw = h * sig[l]
                dens[t, l] += (np.exp(-0.5 * (zd / bw) ** 2) / (np.sqrt(2 * np.pi) * bw)).sum(0)
        done += m
    EGa /= N; Eg /= N; Ea /= N; dens /= N
    gap = float(np.max(np.abs(Eg - gate_ref)))
    print(f"{path}: same-sample check max|P(z>0) - atlas gate_p| = {gap:.2e}", flush=True)
    Gam = EGa - Eg[:, :, None] * Ea[:, None, :]
    np.savez(out, Gam=Gam, w2hat=dens, hs=np.array(hs), Eg=Eg, Ea=Ea, gate_gap=gap)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], k56="--k56" in sys.argv)
