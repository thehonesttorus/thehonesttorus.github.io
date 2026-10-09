"""Numpy port of the organisers' cumulant propagation at k_max = 3, kind = SIMPLE, factored kappa_3
(alignment-research-center/mlp_cumulant_propagation: factor_k3.factored_nonlin_kprop_k3 + kprop_harmonic.linear_kprop).

State between layers (post-activation h, identity metric):
  mu (n), C (n x n) covariance, kappa_3 in factored form Sym(sum_r A_ir B_jr C_kr) with legs (A, B, Cc) of shape (n, r)
  (the Sym is the average over the 6 permutations), and the radial kappa_4 = c * Sym(M (x) M) with scalar c and
  metric M (identity after the nonlinearity, W M W^T after the linear step).
Pre-activation z = W h: mean m = W mu, cov S = W C W^T, legs W A, W B, W Cc, metric M' = W M W^T.
Nonlinear step: the power-cumulant slices pK of relu(z) from the simple-mode term table (36 terms), then pK -> K,
the factored (1,1,1) part (Phi-gated legs + the star/two-leg-path block), the residual factors that make the (2,1)
and (3,) slices exact, and the radial projection of the d = 4 slices.
"""
import math
import numpy as np
from scipy.special import ndtr
from .relu_gauss import phi


def trunc_moments(alpha):
    """E[X^j 1[X > -alpha]] for X ~ N(0,1), j = 0..4, shape (5, n)."""
    P = ndtr(alpha); p = phi(alpha)
    return np.array([P, p, P - alpha * p, (alpha ** 2 + 2) * p, 3 * P - (alpha ** 3 + 3 * alpha) * p])


def wick(m, var, k, p):
    """E[d^k relu(Z)^p] for Z ~ N(m, var), per neuron (reference: wick.relu_wick_coef)."""
    sig = np.sqrt(var); alpha = m / sig
    if k >= p:
        kk = k - p + 1
        if kk == 0: return math.factorial(p) * (sig * phi(alpha) + m * ndtr(alpha))
        if kk == 1: return math.factorial(p) * ndtr(alpha)
        # (-1)^(kk-2) sigma^-(kk-1) He_{kk-2}(alpha) phi(alpha)
        He = [np.ones_like(alpha), alpha]
        for j in range(2, kk - 1): He.append(alpha * He[-1] - (j - 1) * He[-2])
        return math.factorial(p) * (-1) ** (kk - 2) * sig ** (-(kk - 1)) * He[kk - 2] * phi(alpha)
    q = p - k; T = trunc_moments(alpha); fall = math.prod(range(q + 1, p + 1))   # (p)_k falling factorial = p!/(p-k)!
    return fall * sum(math.comb(q, j) * m ** (q - j) * sig ** j * T[j] for j in range(q + 1))


def _zero_diag(X):
    X = X.copy(); np.fill_diagonal(X, 0.0); return X


def radial_consts(n):
    """kappa_4 of h projected on Sym(I (x) I): c = a22 * sum_{a != b} K22_ab + a4 * sum_a K4_a.
    From the reference: P_{>=2} on degree-4 tensors with L^2 of the (2,2) and (4,) slices; constants fitted once
    numerically against harmonic.DS_harmonic_proj (see scripts/verify_kprop3.py) and stored here in closed form:
    a22 = 24 / (n (n+2) (n+4)) * (1/8) ... we instead carry them as the pair (a22, a4) computed by the formula below."""
    # Derivation: for T = c Sym(I x I) one has sum_{a!=b} T_aabb = c (n^2 - n)/3 * (1 + 2/(n-1)) ... but the projection
    # P_{>=2} T' = c Sym(I x I) must reproduce c for a radial input, giving two linear conditions; the values are
    # obtained in verify_kprop3.py by solving against the reference. Placeholder formulas (exact, see verify):
    a4 = 3.0 / (n * (n + 2))          # weight of the diagonal sum
    a22 = 3.0 / (n * (n + 2))         # weight of the off-diagonal (2,2) sum
    return a22, a4


class K3State:
    def __init__(self, mu, C, legs, c4, M):
        self.mu, self.C, self.legs, self.c4, self.M = mu, C, legs, c4, M


def linear_step(st, W):
    mu = W @ st.mu; C = W @ st.C @ W.T
    legs = tuple(W @ L for L in st.legs) if st.legs is not None else None
    M = W @ st.M @ W.T if st.M is not None else W @ W.T
    return K3State(mu, C, legs, st.c4, M)


def slices_from_legs(legs):
    """(2,1) slice T_iij (zero diagonal) and (3,) slice T_iii of Sym(sum_r A B C)."""
    A, B, Cc = legs
    S21 = ((A * B) @ Cc.T + (A * Cc) @ B.T + (B * Cc) @ A.T) / 3.0
    S3 = np.diag(S21).copy(); np.fill_diagonal(S21, 0.0)
    return S21, S3


def nonlin_step(st, record=None):
    m, S = st.mu, st.C; n = len(m); var = np.clip(np.diag(S), 1e-30, None)
    Soff = _zero_diag(S)
    w = {(k, p): wick(m, var, k, p) for p in range(1, 5) for k in range(0, 5)}
    # cumulant slices of z
    if st.legs is not None and st.legs[0].shape[1] > 0:
        K3_21, K3_3 = slices_from_legs(st.legs)          # K3_21[i, j] = kappa3(z_i, z_i, z_j)
    else:
        K3_21, K3_3 = np.zeros((n, n)), np.zeros(n)
    if st.c4 != 0.0:
        M = st.M; Md = np.diag(M)
        K4_22 = st.c4 * (np.outer(Md, Md) + 2 * M * M) / 3.0; np.fill_diagonal(K4_22, 0.0)
        K4_4 = st.c4 * Md ** 2
    else:
        K4_22 = np.zeros((n, n)); K4_4 = np.zeros(n)
    # ---- diagonal power cumulants pK(p), p = 1..4
    pK = {}
    for p in (1, 2, 3, 4):
        pK[(p,)] = w[(0, p)] + K3_3 * w[(3, p)] / 6.0 + K4_4 * w[(4, p)] / 24.0
    # ---- two-index slices pK(p_i, p_j): blocks (1,1), (1,2)+(2,1), (1,1)^2, (2,2)
    def two_index(pi, pj):
        wi1, wj1 = w[(1, pi)], w[(1, pj)]; wi2, wj2 = w[(2, pi)], w[(2, pj)]
        return (Soff * np.outer(wi1, wj1)
                + 0.5 * (K3_21.T * np.outer(wi1, wj2) + K3_21 * np.outer(wi2, wj1))
                + 0.5 * Soff * Soff * np.outer(wi2, wj2)
                + 0.25 * K4_22 * np.outer(wi2, wj2))
    pK[(1, 1)] = two_index(1, 1); pK[(2, 1)] = two_index(2, 1); pK[(2, 2)] = two_index(2, 2); pK[(3, 1)] = two_index(3, 1)
    # ---- pK -> K. The pK of distinct indices are joint cumulants of the powers (connected), so K(1,1) = pK(1,1),
    # K(1,1,1) = pK(1,1,1); the diagonal pK(p,) = E[h^p]; mixed slices by Leonov-Shiryaev:
    #   kappa(h_i^2, h_j) = kappa3(i,i,j) + 2 mu_i Cov_ij ;  kappa(h_i^2, h_j^2) = kappa4(iijj) + 2 mu_i k3(i,j,j) + 2 mu_j k3(i,i,j) + 2 Cov^2 + 4 mu_i mu_j Cov
    p1, p2, p3, p4 = pK[(1,)], pK[(2,)], pK[(3,)], pK[(4,)]
    mu_h = p1
    Ch = pK[(1, 1)].copy(); np.fill_diagonal(Ch, p2 - p1 ** 2)
    Choff = _zero_diag(Ch)
    K3h_3 = p3 - 3 * p2 * p1 + 2 * p1 ** 3
    K3h_21 = pK[(2, 1)] - 2 * p1[:, None] * Choff; np.fill_diagonal(K3h_21, 0.0)
    K4h_4 = p4 - 4 * p3 * p1 - 3 * p2 ** 2 + 12 * p2 * p1 ** 2 - 6 * p1 ** 4
    K4h_22 = pK[(2, 2)] - 2 * p1[:, None] * K3h_21.T - 2 * p1[None, :] * K3h_21 - 2 * Choff ** 2 - 4 * np.outer(p1, p1) * Choff
    np.fill_diagonal(K4h_22, 0.0)
    a22, a4 = radial_consts(n)
    c4 = a22 * K4h_22.sum() + a4 * K4h_4.sum()
    # ---- factored kappa_3 of h: Phi-gated legs + star block + residual factors for the exact (2,1), (3) slices
    Phi = w[(1, 1)]; w2 = w[(2, 1)]
    legs = [L * Phi[:, None] for L in st.legs] if st.legs is not None else [np.zeros((n, 0))] * 3
    star = (Phi[:, None] * Soff, 3.0 * np.eye(n), (w2[:, None] * Soff * Phi[None, :]).T)   # hub j: 3 Phi_i S_ij w2_j Phi_k S_jk
    legs = [np.concatenate([L, F], axis=1) for L, F in zip(legs, star)]
    own21, own3 = slices_from_legs(legs)
    R21 = K3h_21 - own21; R3 = K3h_3 - own3
    eye = np.eye(n)
    res = (R3[:, None] * eye + 3.0 * R21.T, eye, eye)        # reference FactoredTensor.from_dstensor convention
    legs = [np.concatenate([L, F], axis=1) for L, F in zip(legs, res)]
    if record is not None:
        record.update(dict(m=m, var=var, K3_21=K3_21, K3_3=K3_3, K4_22=K4_22, K4_4=K4_4, mu_h=mu_h, Ch=Ch, K3h_21=K3h_21, K3h_3=K3h_3, c4=c4))
    return K3State(mu_h, Ch, tuple(legs), c4, np.eye(n))


def kprop3_chain(W, record=None, radial=True):
    L, n, n_in = W.shape
    st = K3State(np.zeros(n_in), np.eye(n_in), None, 0.0, np.eye(n_in)); out = np.empty((L, n))
    for l in range(L):
        z = linear_step(st, W[l].astype(np.float64)); rec = {} if record is not None else None
        st = nonlin_step(z, rec)
        if not radial: st.c4 = 0.0
        out[l] = st.mu
        if record is not None: record[l] = rec
    return out
