# Is the pre-activation fourth cumulant one coherent product core? (note XXXI, THEORY_3 sections 3 and 7)
#   python k4coherence.py NET MCPREFIX        (MCPREFIX_{full,h0,h1}.npz from mcstats.py, chain_off{NET}.npz from oracle_one dump)
# THEORY_3: a declared core J(M, N) (polynomial (x^T M x)(x^T N x)) fixes every slice at once,
#   d_i = M_ii N_ii,   B_ij = kappa(i,i,i,j) = (M_ii N_ij + M_ij N_ii)/2,   K_ij = kappa(i,i,j,j) = (M_ii N_jj + M_jj N_ii + 4 M_ij N_ij)/6.
# Given a metric M, the diagonal and the (3,1) slice determine N (N_ii = d_i / M_ii, N_ij = (2 B_ij - M_ij N_ii) / M_ii);
# the (2,2) slice is then a PREDICTION, and N_ij against N_ji tests the (3,1) orientation. Metrics tested:
#   C   = the true covariance of z_l (the scale mixture is J(C, 3gC)); exactly covariant under the linear layer;
#   WW  = W_l W_l^T (THEORY_3's transported Euclidean trace metric of the previous post-activation);
#   2I  = the chain's METRIC_C convention.
# Applied to Monte Carlo truth (is the physical kappa_4 one core, and in which metric?) and to the chain's own declared slices
# (g4row, wk4m, wk431: are they one tensor?). Off-diagonal statistics only; MC noise from the two halves.
import sys, numpy as np
net, pre = int(sys.argv[1]), sys.argv[2]
F, H0, H1 = (np.load(f"{pre}_{t}.npz") for t in ("full", "h0", "h1"))
ch = np.load(f"chain_off{net}.npz")
Wcol = np.load(f"../official/W_off{net}.npy").astype(np.float64)
n = Wcol.shape[1]; off = ~np.eye(n, dtype=bool)


def rel(p, t):
    return np.linalg.norm(p - t) / np.linalg.norm(t)


def core_predict(M, d, B):
    """N from diagonal and (3,1) slice in metric M; returns predicted K22 (off-diag) and the asymmetry of N."""
    m = np.diag(M)
    N = (2.0 * B - M * (d / m)[:, None]) / m[:, None]          # row i: from B_ij = kappa(i,i,i,j)
    np.fill_diagonal(N, d / m)
    asym = rel(N[off], N.T[off]) if np.linalg.norm(N[off]) > 0 else 0.0
    Ns = 0.5 * (N + N.T); nd = np.diag(Ns)
    K = (m[:, None] * nd[None, :] + m[None, :] * nd[:, None] + 4.0 * M * Ns) / 6.0
    return K, asym


print(f"net {net}: {int(F['n'])} samples. Per layer: MC noise of the slices, then for each metric the (2,2) slice predicted from the"
      f" diagonal + (3,1) slice of one core J(M, N) (rel err, corr vs MC K22) and the asymmetry of N; mixture J(C, 3gC) for reference")
for l in range(2, 16):
    C = F["cov"][l].astype(np.float64); var = F["var"][l].astype(np.float64); np.fill_diagonal(C, var)
    d = F["k4"][l].astype(np.float64); K = F["K22"][l].astype(np.float64); B = F["K31"][l].astype(np.float64)
    nzK = np.linalg.norm((H0["K22"][l] - H1["K22"][l])[off]) / 2 / np.linalg.norm(K[off])
    nzB = np.linalg.norm((H0["K31"][l] - H1["K31"][l])[off]) / 2 / np.linalg.norm(B[off])
    nzd = np.linalg.norm(H0["k4"][l] - H1["k4"][l]) / 2 / np.linalg.norm(d)
    W = Wcol[l]
    mets = {"C": C, "WW": W @ W.T, "2I": 2.0 * np.eye(n)}
    parts = []
    for nm, M in mets.items():
        Kp, asym = core_predict(M, d, B)
        parts.append(f"{nm}: K22 rel {rel(Kp[off], K[off]):.2f} corr {np.corrcoef(Kp[off], K[off])[0, 1]:+.3f} asym {asym:.2f}")
    g = float(np.sum(d * var ** 2) / np.sum(3 * var ** 4))
    Kmix = g * (var[:, None] * var[None, :] + 2 * C * C); Bmix = 3 * g * var[:, None] * C
    mix = (f"mixture g={g:.4f}: diag rel {rel(3 * g * var ** 2, d):.2f} | K22 rel {rel(Kmix[off], K[off]):.2f} corr "
           f"{np.corrcoef(Kmix[off], K[off])[0, 1]:+.3f} | K31 rel {rel(Bmix[off], B[off]):.2f} corr {np.corrcoef(Bmix[off], B[off])[0, 1]:+.3f}")
    # the chain's own declared slices: one tensor?
    cd = ch[f"g4row_{l}"].astype(np.float64); cK = ch[f"wk4m_{l}"].astype(np.float64)
    cB = ch[f"wk431_{l}"].astype(np.float64).T                  # wk431[a, c] = kappa(a, c, c, c)
    cC = ch[f"C_off_{l}"].astype(np.float64).copy(); cv = ch[f"var_{l}"].astype(np.float64); np.fill_diagonal(cC, cv)
    Kc, asc = core_predict(cC, cd, cB)
    chain = (f"chain's own slices in metric C: predicted K22 vs its own wk4m rel {rel(Kc[off], cK[off]):.2f} corr "
             f"{np.corrcoef(Kc[off], cK[off])[0, 1]:+.3f} asym {asc:.2f}; chain vs MC: diag {rel(cd, d):.2f} K22 {rel(cK[off], K[off]):.2f} "
             f"K31 {rel(cB[off], B[off]):.2f}")
    print(f"layer {l:2d} (MC noise diag {nzd:.3f} K22 {nzK:.3f} K31 {nzB:.3f}):\n    " + " | ".join(parts) + f"\n    {mix}\n    {chain}",
          flush=True)
