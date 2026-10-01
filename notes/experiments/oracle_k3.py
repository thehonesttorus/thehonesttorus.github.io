"""Oracle tests of third-cumulant representations on a small-width atlas (numpy only).

Input: atlas npz files from moment_atlas_np.py --k3 (full raw third-moment tensors pre_M3, post_M3).
For each layer l it forms the exact cumulant objects the K=3 chains consume and tests how well
compressed representations of the post-activation third cumulant kappa3(a_l) reproduce the next
layer's pre-activation (2,1) slice D21(l+1)_{ab} = kappa3(z_a, z_a, z_b), which the published chain
established as the whole interface to the next nonlinearity (error law: extra MSE ~ 4.2e-6 eps^2,
eps = relative rms error of D21; eps <= 2.2 % needed at the frontier).

Single-atlas table (python oracle_k3.py A.npz [B.npz ...]), one row per layer:
  transport   relative error of D21(l+1) computed by contracting kappa3(a_l) with W_{l+1} against the
              atlas's own D21(l+1): an identity at the sample level, so it only checks the code
  R2(K22|C)   R^2 of the off-diagonal (2,2) pre-activation cumulant slice on C_off, and on C_off + C_off^2
  d21 n/8..   fraction of ||D21||_F^2 captured by rank n/8, n/4, n/2 (SVD)
  memless     eps when kappa3(a) is replaced by its (3,) and (2,1) slices only (504aldo's dead end 2);
              'off' = the same eps restricted to the a != b entries of D21
  wick        eps of the leading-order Wick closure of the all-distinct part (the term every factorised
              K=3 chain carries implicitly in O(n^3)):
                 kappa3(a)_{ijk} ~ sum_cyc Phi_i Phi_j w2_k C_ik C_jk + Phi_i Phi_j Phi_k kappa3(z)_{ijk}
              with Phi = P(z > 0) (atlas gate_p) and w2 = phi(alpha)/sigma (Gaussian density of z at 0)
  hub n/r     eps when the all-distinct part is kept as the top n/r "hub columns" (SVD of the (n^2, n) unfolding)
  sketch s    eps when the all-distinct part is projected on a random s-dimensional subspace along two indices
  wick+hub    eps when the Wick closure is used and its residual is kept as top n/4 hub columns

Pair mode (python oracle_k3.py --pair A.npz B.npz) needs two atlases of the SAME MLP built with
different --sample-seed; it measures the Monte Carlo noise floor of every quantity and re-evaluates
the representations fitted on A against the independent target of B, noise-corrected:
  eps_noise   relative noise of one atlas's D21 (all entries / off-diagonal), = rel(D21_A, D21_B)/sqrt(2)
  snr         signal-to-noise energy ratio of the all-distinct kappa3(a) tensor, <A,B>/(mean(|A|^2,|B|^2) - <A,B>)
  memless x   eps of A's slices-only model against B's D21, raw and after subtracting B's noise energy
  wick x      the same for the Wick closure built from A
  hub x       the same for A's top-n/4 hub columns (subspace and coefficients from A)
"""
import sys
import numpy as np
from math import erf, sqrt, pi, exp


def central3(M3, M2, m1):
    """kappa3 tensor from raw moments: E[(x-mu)_i (x-mu)_j (x-mu)_k]."""
    t = M3.copy()
    t -= np.einsum("i,jk->ijk", m1, M2) + np.einsum("j,ik->ijk", m1, M2) + np.einsum("k,ij->ijk", m1, M2)
    t += 2.0 * np.einsum("i,j,k->ijk", m1, m1, m1)
    return t


def slices(K3):
    n = K3.shape[0]
    idx = np.arange(n)
    return K3[idx, idx, idx], K3[idx, idx, :]      # D3[a] = kappa3(a,a,a), D21[a, b] = kappa3(a,a,b)


def rel(a, b, off=False):
    if off:
        a, b = offdiag(a), offdiag(b)
    return float(np.sqrt(np.sum((a - b) ** 2) / np.sum(b ** 2)))


def offdiag(M):
    M = M.copy(); np.fill_diagonal(M, 0.0); return M


def r2_fit(y, X):
    """uncentred R^2 of least squares y ~ X (columns), both flattened."""
    X = np.stack([x.ravel() for x in X], 1); yv = y.ravel()
    coef, *_ = np.linalg.lstsq(X, yv, rcond=None)
    res = yv - X @ coef
    return 1.0 - float(res @ res) / float(yv @ yv), coef


def transport_d21(K3a, W):
    """D21(l+1)_{ab} = sum_{ijk} W_ia W_ja W_kb kappa3(a_i,a_j,a_k)."""
    T = np.einsum("ijk,kb->ijb", K3a, W)
    T = np.einsum("ijb,ja->iab", T, W)
    return np.einsum("iab,ia->ab", T, W)


def all_distinct(K3):
    n = K3.shape[0]
    K = K3.copy()
    idx = np.arange(n)
    K[idx, idx, :] = 0.0; K[idx, :, idx] = 0.0; K[:, idx, idx] = 0.0
    return K


def slices_only(K3):
    """the (3,) and (2,1) slices of K3 in tensor form (zero on all-distinct entries)."""
    return K3 - all_distinct(K3)


def hub_columns(Kd, k):
    n = Kd.shape[0]
    U, S, Vt = np.linalg.svd(Kd.reshape(n * n, n), full_matrices=False)
    return ((U[:, :k] * S[:k]) @ Vt[:k]).reshape(n, n, n)


def wick_model(C, Phi, w2, K3z):
    """Leading-order Wick closure of the all-distinct post-activation third cumulant (see module docstring)."""
    T = np.einsum("i,j,k,ik,jk->ijk", Phi, Phi, w2, C, C)
    T += np.einsum("i,j,k,ij,kj->ijk", Phi, w2, Phi, C, C)
    T += np.einsum("i,j,k,ji,ki->ijk", w2, Phi, Phi, C, C)
    T += np.einsum("i,j,k,ijk->ijk", Phi, Phi, Phi, K3z)
    return all_distinct(T)


def gaussian_w2(mu, var):
    sig = np.sqrt(var); alpha = mu / sig
    return np.exp(-0.5 * alpha ** 2) / np.sqrt(2 * pi) / sig


def layer_objects(z, l):
    """the cumulant objects of layer l (pre-activation z_l and post-activation a_l) and the next layer's D21."""
    pre_s, post_s = z["pre_s"], z["post_s"]
    mu, m2 = pre_s[0, l], pre_s[1, l]; var = m2 - mu ** 2
    M11 = z["pre_M11"][l].astype(np.float64)
    C = M11 - np.outer(mu, mu)
    K3z = central3(z["pre_M3"][l], M11, mu)
    mu_a = post_s[0, l]
    K3a = central3(z["post_M3"][l], z["post_M11"][l].astype(np.float64), mu_a)
    mu1 = pre_s[0, l + 1]
    K3z1 = central3(z["pre_M3"][l + 1], z["pre_M11"][l + 1].astype(np.float64), mu1)
    D3_1, D21_1 = slices(K3z1)
    Phi = z["gate_p"][l].astype(np.float64)
    return dict(mu=mu, var=var, C=C, K3z=K3z, K3a=K3a, D21=D21_1, Phi=Phi, w2=gaussian_w2(mu, var))


def k22_lambda_law(z, l1):
    """R^2 of the off-diagonal (2,2) cumulant slice of the pre-activation at layer l1 on C_off and C_off^2."""
    pre_s = z["pre_s"]
    mu1 = pre_s[0, l1]; m2_1 = pre_s[1, l1]; var1 = m2_1 - mu1 ** 2
    M11 = z["pre_M11"][l1].astype(np.float64); M21 = z["pre_M21"][l1].astype(np.float64); M22 = z["pre_M22"][l1].astype(np.float64)
    C1 = M11 - np.outer(mu1, mu1)
    Eu2u2 = (M22 - 2 * mu1[None, :] * M21 - 2 * mu1[:, None] * M21.T
             + np.outer(m2_1, mu1 ** 2) + np.outer(mu1 ** 2, m2_1) + 4 * np.outer(mu1, mu1) * M11
             - 3 * np.outer(mu1 ** 2, mu1 ** 2))
    K22 = Eu2u2 - np.outer(var1, var1) - 2 * C1 ** 2
    K22o, Co = offdiag(K22), offdiag(C1)
    return r2_fit(K22o, [Co])[0], r2_fit(K22o, [Co, Co * Co])[0]


def analyse(path, ranks=(8, 4, 2), sketch=(8, 16, 32), rng=np.random.default_rng(0)):
    z = np.load(path)
    W = z["weights"].astype(np.float64)
    L, n, _ = W.shape
    print(f"\n{path}: width {n}, depth {L}, N = {int(z['n_samples'])}")
    print(f"{'l':>2} {'transp':>6} {'R2(K22|C)':>9} {'R2+C*C':>7} {'d21n/8':>6} {'n/4':>5} {'n/2':>5} | {'memless':>7} {'off':>5} | {'wick':>5} {'off':>5} | "
          + " ".join(f"hub{n//r:>3}" for r in ranks) + " | " + " ".join(f"sk{s:>3}" for s in sketch) + f" | {'wick+hub':>8}")
    for l in range(L - 1):
        o = layer_objects(z, l)
        K3a, D21 = o["K3a"], o["D21"]
        Wn = W[l + 1]
        e_transport = rel(transport_d21(K3a, Wn), D21)
        r2a, r2b = k22_lambda_law(z, l + 1)
        sv = np.linalg.svd(D21, compute_uv=False); e = sv ** 2 / np.sum(sv ** 2)
        spec = [float(np.sum(e[: n // r])) for r in ranks]
        K3m = slices_only(K3a)
        Tm = transport_d21(K3m, Wn)
        e_mem, e_mem_off = rel(Tm, D21), rel(Tm, D21, off=True)
        Kw = wick_model(o["C"], o["Phi"], o["w2"], o["K3z"])
        Tw = transport_d21(K3m + Kw, Wn)
        e_wick, e_wick_off = rel(Tw, D21), rel(Tw, D21, off=True)
        Kd = all_distinct(K3a)
        hub = [rel(transport_d21(K3m + hub_columns(Kd, n // r), Wn), D21) for r in ranks]
        sk = []
        for s_ in sketch:
            Om = rng.standard_normal((n, s_)) / np.sqrt(s_)
            Q, _ = np.linalg.qr(Om)
            P = Q @ Q.T
            Ks = np.einsum("ijk,jm,kn->imn", Kd, P, P)
            sk.append(rel(transport_d21(K3m + Ks, Wn), D21))
        e_wh = rel(transport_d21(K3m + Kw + hub_columns(Kd - Kw, n // 4), Wn), D21)
        print(f"{l:>2} {e_transport:6.3f} {r2a:9.3f} {r2b:7.3f} {spec[0]:6.3f} {spec[1]:5.3f} {spec[2]:5.3f} | {e_mem:7.3f} {e_mem_off:5.3f} | {e_wick:5.3f} {e_wick_off:5.3f} | "
              + " ".join(f"{h:6.3f}" for h in hub) + " | " + " ".join(f"{v:5.3f}" for v in sk) + f" | {e_wh:8.3f}")


def corrected(e, e_noise):
    return float(np.sqrt(max(e ** 2 - e_noise ** 2, 0.0)))


def analyse_pair(pa, pb):
    A, B = np.load(pa), np.load(pb)
    assert np.array_equal(A["weights"], B["weights"]), "pair mode needs two atlases of the same MLP"
    W = A["weights"].astype(np.float64)
    L, n, _ = W.shape
    print(f"\npair {pa} | {pb}: width {n}, depth {L}, N = {int(A['n_samples'])} + {int(B['n_samples'])}")
    print(f"{'l':>2} {'eps_noise':>9} {'off':>6} {'snr(Kd)':>8} | {'memless x':>9} {'corr':>6} {'off':>6} {'corr':>6} | {'wick x':>7} {'corr':>6} {'off':>6} {'corr':>6} | {'hub n/4 x':>9} {'corr':>6}")
    for l in range(L - 1):
        oa, ob = layer_objects(A, l), layer_objects(B, l)
        Wn = W[l + 1]
        D21a, D21b = oa["D21"], ob["D21"]
        e_noise = rel(D21a, D21b) / sqrt(2)
        e_noise_off = rel(D21a, D21b, off=True) / sqrt(2)
        Kda, Kdb = all_distinct(oa["K3a"]), all_distinct(ob["K3a"])
        cross = float(np.sum(Kda * Kdb)); self_ = 0.5 * (float(np.sum(Kda ** 2)) + float(np.sum(Kdb ** 2)))
        snr = cross / max(self_ - cross, 1e-300)
        K3m = slices_only(oa["K3a"])
        Tm = transport_d21(K3m, Wn)
        e_mem, e_mem_off = rel(Tm, D21b), rel(Tm, D21b, off=True)
        Kw = wick_model(oa["C"], oa["Phi"], oa["w2"], oa["K3z"])
        Tw = transport_d21(K3m + Kw, Wn)
        e_wick, e_wick_off = rel(Tw, D21b), rel(Tw, D21b, off=True)
        e_hub = rel(transport_d21(K3m + hub_columns(Kda, n // 4), Wn), D21b)
        print(f"{l:>2} {e_noise:9.3f} {e_noise_off:6.3f} {snr:8.2f} | {e_mem:9.3f} {corrected(e_mem, e_noise):6.3f} {e_mem_off:6.3f} {corrected(e_mem_off, e_noise_off):6.3f} | "
              f"{e_wick:7.3f} {corrected(e_wick, e_noise):6.3f} {e_wick_off:6.3f} {corrected(e_wick_off, e_noise_off):6.3f} | {e_hub:9.3f} {corrected(e_hub, e_noise):6.3f}")


def selftest(n=3, N=4_000_000, seed=0):
    """check the Wick closure's leading term on a correlated Gaussian triple with small correlations."""
    rng = np.random.default_rng(seed)
    A = rng.standard_normal((n, n)) / np.sqrt(n) * 0.25 + np.eye(n)
    C = A @ A.T
    mu = np.array([0.3, -0.2, 0.1])
    z = rng.multivariate_normal(mu, C, size=N)
    a = np.maximum(z, 0.0)
    ac = a - a.mean(0)
    K3a = np.einsum("mi,mj,mk->ijk", ac, ac, ac) / N
    Phi = (z > 0).mean(0)
    w2 = gaussian_w2(mu, np.diag(C))
    Kw = wick_model(C, Phi, w2, np.zeros((n, n, n)))
    print("selftest: all-distinct kappa3(a) by MC", K3a[0, 1, 2], " Wick leading order", Kw[0, 1, 2],
          f" (relative gap {abs(K3a[0,1,2]-Kw[0,1,2])/abs(K3a[0,1,2]):.3f}; MC se ~ {np.sqrt(np.var(ac[:,0]*ac[:,1]*ac[:,2])/N):.2e})")


if __name__ == "__main__":
    args = sys.argv[1:]
    if args and args[0] == "--selftest":
        selftest()
    elif args and args[0] == "--pair":
        analyse_pair(args[1], args[2])
    else:
        for p in args:
            analyse(p)
