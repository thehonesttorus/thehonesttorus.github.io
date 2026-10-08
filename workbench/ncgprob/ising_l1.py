# The exact layer-1 Ising conditional of the signs given the magnitudes (the "exact Gibbs model inside the setting" of the
# pasted synthesis): Z = W1 x ~ N(0, Sigma), Sigma = W1 W1^T, p(sigma | r) ~ exp(sum_{i<j} J_ij sigma_i sigma_j) with
# J_ij(r) = -r_i (Sigma^-1)_ij r_j.  Sizes of the couplings at typical magnitudes and the Dobrushin row sums.
import numpy as np
for net in (0, 1):
    W = np.load(f"../official/W_off{net}.npy")[0].astype(np.float64); n = W.shape[0]
    S = W @ W.T; ev = np.linalg.eigvalsh(S); Si = np.linalg.inv(S)
    r = np.sqrt(np.diag(S))
    J = -np.outer(r, r) * Si; np.fill_diagonal(J, 0.0)
    rowsum = np.sum(np.tanh(np.abs(J)), axis=1); off = J[~np.eye(n, dtype=bool)]
    print(f"net {net}: Sigma eigenvalues min {ev[0]:.3e} max {ev[-1]:.3e} (condition {ev[-1]/ev[0]:.2e}); diag Sigma mean {np.mean(np.diag(S)):.3f}")
    print(f"   |J_ij| at typical magnitudes: rms {np.sqrt(np.mean(off**2)):.3f}, median {np.median(np.abs(off)):.3f}, 99th pct {np.quantile(np.abs(off), 0.99):.2f}, max {np.abs(off).max():.1f}")
    print(f"   Dobrushin row sums sum_j tanh|J_ij|: min {rowsum.min():.1f}, median {np.median(rowsum):.1f}, max {rowsum.max():.1f}  (uniqueness needs < 1)")
    x = np.random.default_rng(net).standard_normal(n); z = W @ x; rr = np.abs(z); Jx = -np.outer(rr, rr) * Si; np.fill_diagonal(Jx, 0.0)
    print(f"   one random input: Dobrushin row sums median {np.median(np.sum(np.tanh(np.abs(Jx)), 1)):.1f}, max {np.sum(np.tanh(np.abs(Jx)), 1).max():.1f}")
