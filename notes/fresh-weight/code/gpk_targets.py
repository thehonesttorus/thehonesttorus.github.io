# Note XLIV section 9 (active-set version: only neurons with alpha = mu / S > -2.5, the set the chain keeps; D3 is NOT masked in the chain): per-layer truth targets for the gain-package oracle (V62_GPK) of the research estimator.
#   python -I gpk_targets.py MC2DIR NET OUT.npz
# e_v = Rayleigh quotient mu^T C mu / |mu|^4 of the true pre-activation covariance; e_3, e_21, e_22, e_31, e_4 = least-squares
# dilation charge of each slice (k3 diagonal, D21, K22, K31, k4 diagonal) on its package form, with the forms built from the true
# mean, variance and covariance (Monte Carlo, both halves averaged); the same definitions as dil_ledger.py and the estimator hook.
import sys, numpy as np
mcd, net, outp = sys.argv[1], int(sys.argv[2]), sys.argv[3]
F = {}
for h in ("h0", "h1"):
    z = np.load(f"{mcd}/mc2_off{net}_{h}.npz")
    F[h] = {k: z[k] for k in ("mu", "var", "k3", "k4", "cov", "D21", "K22", "K31")}
n = F["h0"]["mu"].shape[1]
off = ~np.eye(n, dtype=bool)
f64 = lambda a: np.asarray(a, dtype=np.float64)
T = {k: np.full(16, np.nan) for k in ("e_v", "e_3", "e_21", "e_22", "e_31", "e_4")}
for s in range(1, 15):
    mu = 0.5 * (f64(F["h0"]["mu"][s]) + f64(F["h1"]["mu"][s])); v = 0.5 * (f64(F["h0"]["var"][s]) + f64(F["h1"]["var"][s]))
    C = 0.5 * (f64(F["h0"]["cov"][s]) + f64(F["h1"]["cov"][s])); C = 0.5 * (C + C.T)
    act = (mu / np.sqrt(v)) > -2.5
    mA = mu * act
    T["e_v"][s] = float(mA @ C @ mA / (mA @ mA) ** 2)
    mk = off & act[:, None] & act[None, :]
    f3 = 6 * mu * v; f4 = 12 * v * v
    f21 = np.where(mk, 2 * (2 * mu[:, None] * C + mu[None, :] * v[:, None]), 0.0)
    f22 = np.where(mk, 4 * np.outer(v, v) + 8 * C * C, 0.0)
    f31 = np.where(mk, 12 * v[:, None] * C, 0.0)            # K31[a, b] = kappa(a, a, a, b) = 12 var_a C_ab
    def ch(k, f, dim):
        if dim == 2:
            return float(np.mean([np.sum(f * f64(F[h][k][s])) / np.sum(f * f) for h in ("h0", "h1")]))
        return float(np.mean([f[act] @ f64(F[h][k][s])[act] / (f[act] @ f[act]) for h in ("h0", "h1")]))
    T["e_3"][s] = ch("k3", f3, 1); T["e_4"][s] = ch("k4", f4, 1)
    T["e_21"][s] = ch("D21", f21, 2); T["e_22"][s] = ch("K22", f22, 2); T["e_31"][s] = ch("K31", f31, 2)
np.savez(outp, **T)
print(f"net {net}: saved {outp}")
for s in range(3, 15):
    print(f"  {s:2d} " + " ".join(f"{k} {1e3 * T[k][s]:6.2f}" for k in T))
