# Oracle inputs that price representations of the (3,1) slice: copies of mc2_off{NET}_{full,h0,h1}.npz with K31
# replaced at every layer by
#   colK : K31 U U^T, its projection on the top-K eigenvectors of the (full-sample) covariance on the column index
#          (what any latent / column-collective model can at most supply), for each K in KLIST;
#   gain : its scale-mixture fit 3 g var_i C_ij (g fitted per layer on the full sample), the gain shape with the
#          var factor the production slice lam C_off lacks.
#   python k31proj.py NET KLIST     e.g.  python k31proj.py 0 16,64
# Writes mc2col{K}_off{NET}_{tag}.npz and mc2gain_off{NET}_{tag}.npz; run oracle_one.py NET K31 FILE TAG on each.
import sys, numpy as np
net = int(sys.argv[1]); KS = [int(x) for x in sys.argv[2].split(",")]
F = {t: np.load(f"mc2_off{net}_{t}.npz") for t in ("full", "h0", "h1")}
D = {t: {k: F[t][k] for k in F[t].files} for t in F}
L = D["full"]["K31"].shape[0]
U = {}; G = {}
for l in range(L):
    C = np.asarray(D["full"]["cov"][l], np.float64)
    ev, V = np.linalg.eigh(C); U[l] = V[:, ::-1]
    var = np.asarray(D["full"]["var"][l], np.float64); S = 3 * var[:, None] * C; np.fill_diagonal(S, 0.0)
    X = np.asarray(D["full"]["K31"][l], np.float64).copy(); np.fill_diagonal(X, 0.0)
    G[l] = float(np.sum(X * S)) / max(float(np.sum(S * S)), 1e-300)
for t in D:
    for K in KS:
        out = dict(D[t]); k31 = np.array(D[t]["K31"], np.float32, copy=True)
        for l in range(L):
            Uk = U[l][:, :K]; X = np.asarray(k31[l], np.float64)
            k31[l] = (X @ Uk @ Uk.T).astype(np.float32)
        out["K31"] = k31
        np.savez(f"mc2col{K}_off{net}_{t}.npz", **out)
    out = dict(D[t]); k31 = np.array(D[t]["K31"], np.float32, copy=True)
    for l in range(L):
        var = np.asarray(D[t]["var"][l], np.float64); C = np.asarray(D[t]["cov"][l], np.float64)
        k31[l] = (G[l] * 3 * var[:, None] * C).astype(np.float32)
    out["K31"] = k31
    np.savez(f"mc2gain_off{net}_{t}.npz", **out)
print("gain per layer:", " ".join(f"{G[l]:.4f}" for l in range(L)))
