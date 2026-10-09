"""Two-tier version of the K=3 simple chain (whest/kprop3.py): young kappa_3 sources (age <= window) are kept as dense
factored legs; older sources are merged into a Tucker tier Sym((Q x Q x Q) S) with a shared orthonormal basis Q (n x k)
and a symmetric core S (k x k x k). The basis is the dominant left subspace of the gated transport cocycle, which is
what the sources of age >= window live in (scripts/diag_cocycle.py). Options: window (ages kept dense), k (basis size).
"""
import numpy as np
from .kprop3 import wick, _zero_diag, radial_consts, slices_from_legs, linear_step as _lin_unused  # noqa


class State:
    def __init__(self, mu, C, young, Q, S, c4, M):
        self.mu, self.C, self.young, self.Q, self.S, self.c4, self.M = mu, C, young, Q, S, c4, M   # young: list of (legs, age)


def tucker_slices(Q, S):
    """(2,1) slice (zero diagonal) and (3,) slice of Sym((Q x Q x Q) S) for a symmetric core S."""
    n, k = Q.shape
    M = np.empty((n, k))
    for s in range(k):
        M[:, s] = np.einsum("ip,ip->i", Q @ S[:, :, s], Q)      # sum_pq Q_ip Q_iq S_pqs
    S21 = M @ Q.T; S3 = np.einsum("is,is->i", M, Q); np.fill_diagonal(S21, 0.0)
    return S21, S3


def sym3(S):
    return (S + S.transpose(0, 2, 1) + S.transpose(1, 0, 2) + S.transpose(1, 2, 0) + S.transpose(2, 0, 1) + S.transpose(2, 1, 0)) / 6.0


def linear_step(st, W):
    mu = W @ st.mu; C = W @ st.C @ W.T
    young = [(tuple(W @ L for L in legs), age) for legs, age in st.young]
    Q, S = st.Q, st.S
    if Q is not None:
        Y = W @ Q; Q, R = np.linalg.qr(Y)
        S = np.einsum("pqs,ap,bq,cs->abc", S, R, R, R, optimize=True)
    return State(mu, C, young, Q, S, st.c4, W @ st.M @ W.T)


def merge_into_tucker(Q, S, legs, k, rng, fixed=None):
    """Add the factored block Sym(sum_r A B C) to the Tucker tier, refitting the basis to the top-k directions of
    [Q diag(nu), A, B, C] with nu the core's mode norms (randomized range finder, one power iteration). `fixed`
    (n x f) directions are always kept in the basis (the collective directions), the rest fills up to k."""
    A, B, Cc = legs; n = A.shape[0]
    if Q is None:
        G = np.concatenate([A, B, Cc], axis=1)
    else:
        nu = np.sqrt(np.einsum("pqs,pqs->p", S, S)); G = np.concatenate([Q * nu[None, :], A, B, Cc], axis=1)
    Om = rng.standard_normal((G.shape[1], k + 8)); Y = G @ Om; Y = G @ (G.T @ Y)
    if fixed is not None and fixed.shape[1] > 0:
        F = np.linalg.qr(fixed)[0]; Y = Y - F @ (F.T @ Y)
        Qn = np.concatenate([F, np.linalg.qr(Y)[0][:, :k - F.shape[1]]], axis=1)
    else:
        Qn = np.linalg.qr(Y)[0][:, :k]
    Sn = np.zeros((k, k, k))
    if Q is not None:
        P = Qn.T @ Q; Sn += np.einsum("pqs,ap,bq,cs->abc", S, P, P, P, optimize=True)
    a, b, c = Qn.T @ A, Qn.T @ B, Qn.T @ Cc
    Sn += np.einsum("ar,br,cr->abc", a, b, c, optimize=True)
    return Qn, sym3(Sn)


def nonlin_step(st, o, rng, record=None):
    m, S_ = st.mu, st.C; n = len(m); var = np.clip(np.diag(S_), 1e-30, None); Soff = _zero_diag(S_)
    w = {(k, p): wick(m, var, k, p) for p in range(1, 5) for k in range(0, 5)}
    K3_21 = np.zeros((n, n)); K3_3 = np.zeros(n); mode = o.get("oldmode", "full")
    for legs, age in st.young:
        s21, s3 = slices_from_legs(legs)
        if age >= o["window"] and mode != "full":          # diagnostic: what the old sources' slices are worth
            if mode == "none": continue
            if mode == "diag": s21 = 0.0 * s21
            if mode.startswith("scale"): f = float(mode[5:]); s21 = f * s21; s3 = f * s3
        K3_21 += s21; K3_3 += s3
    if st.Q is not None:
        s21, s3 = tucker_slices(st.Q, st.S); K3_21 += s21; K3_3 += s3
    if st.c4 != 0.0:
        Md = np.diag(st.M); K4_22 = st.c4 * (np.outer(Md, Md) + 2 * st.M * st.M) / 3.0; np.fill_diagonal(K4_22, 0.0); K4_4 = st.c4 * Md ** 2
    else:
        K4_22 = np.zeros((n, n)); K4_4 = np.zeros(n)
    pK = {}
    for p in (1, 2, 3, 4):
        pK[(p,)] = w[(0, p)] + K3_3 * w[(3, p)] / 6.0 + K4_4 * w[(4, p)] / 24.0
    def two_index(pi, pj):
        wi1, wj1 = w[(1, pi)], w[(1, pj)]; wi2, wj2 = w[(2, pi)], w[(2, pj)]
        return (Soff * np.outer(wi1, wj1) + 0.5 * (K3_21.T * np.outer(wi1, wj2) + K3_21 * np.outer(wi2, wj1))
                + 0.5 * Soff * Soff * np.outer(wi2, wj2) + 0.25 * K4_22 * np.outer(wi2, wj2))
    pK[(1, 1)] = two_index(1, 1); pK[(2, 1)] = two_index(2, 1); pK[(2, 2)] = two_index(2, 2)
    p1, p2, p3, p4 = pK[(1,)], pK[(2,)], pK[(3,)], pK[(4,)]
    mu_h = p1; Ch = pK[(1, 1)].copy(); np.fill_diagonal(Ch, p2 - p1 ** 2); Choff = _zero_diag(Ch)
    K3h_3 = p3 - 3 * p2 * p1 + 2 * p1 ** 3
    K3h_21 = pK[(2, 1)] - 2 * p1[:, None] * Choff; np.fill_diagonal(K3h_21, 0.0)
    K4h_4 = p4 - 4 * p3 * p1 - 3 * p2 ** 2 + 12 * p2 * p1 ** 2 - 6 * p1 ** 4
    K4h_22 = pK[(2, 2)] - 2 * p1[:, None] * K3h_21.T - 2 * p1[None, :] * K3h_21 - 2 * Choff ** 2 - 4 * np.outer(p1, p1) * Choff
    np.fill_diagonal(K4h_22, 0.0)
    a22, a4 = radial_consts(n); c4 = (a22 * K4h_22.sum() + a4 * K4h_4.sum()) * o.get("c4scale", 1.0)
    Phi = w[(1, 1)]; w2 = w[(2, 1)]
    # gate all tiers, build the star block, then the residual block from the exact slices
    young = [(tuple(L * Phi[:, None] for L in legs), age + 1) for legs, age in st.young]
    Q, S = st.Q, st.S
    if Q is not None:
        Qg, R = np.linalg.qr(Phi[:, None] * Q); Q = Qg; S = np.einsum("pqs,ap,bq,cs->abc", S, R, R, R, optimize=True)
    star = (Phi[:, None] * Soff, 3.0 * np.eye(n), (w2[:, None] * Soff * Phi[None, :]).T)
    # slices of the gated transported sources = gates applied to the slices already read (Phi_i^2 Phi_j, Phi_i^3)
    own21 = (Phi ** 2)[:, None] * K3_21 * Phi[None, :]; own3 = Phi ** 3 * K3_3
    s21, s3 = slices_from_legs(star); own21 = own21 + s21; own3 += s3
    R21 = K3h_21 - own21; R3 = K3h_3 - own3; eye = np.eye(n)
    newborn = tuple(np.concatenate([F, G], axis=1) for F, G in zip(star, (R3[:, None] * eye + 3.0 * R21.T, eye, eye)))
    young.append((newborn, 0))
    # aging: merge sources of age >= window into the Tucker tier
    keep = []; fixed = None
    if o.get("collective", 0):
        v = mu_h.copy()
        for _ in range(4): v = Ch @ v; v /= np.linalg.norm(v)          # top eigenvector of the post-activation covariance
        fixed = np.stack([mu_h, np.diag(Ch), np.ones(n), v], axis=1)
    for legs, age in young:
        if age >= o["window"] and o.get("tucker", 1):
            Q, S = merge_into_tucker(Q, S, legs, o["k"], rng, fixed)
        else:
            keep.append((legs, age))
    if record is not None:
        record.update(dict(m=m, var=var, K3_21=K3_21, K3_3=K3_3, mu_h=mu_h, Ch=Ch, K3h_21=K3h_21, K3h_3=K3h_3, c4=c4, young=len(keep)))
    return State(mu_h, Ch, keep, Q, S, c4, np.eye(n))


def kprop3c_chain(W, opts=None, record=None, m0=None, S0=None):
    o = dict(window=2, k=128, c4scale=1.0, seed=0, tucker=1, oldmode="full"); o.update(opts or {})
    L, n, n_in = W.shape; rng = np.random.default_rng(o["seed"])
    m0 = np.zeros(n_in) if m0 is None else np.asarray(m0, dtype=np.float64); S0 = np.eye(n_in) if S0 is None else np.asarray(S0, dtype=np.float64)
    st = State(m0.copy(), S0.copy(), [], None, None, 0.0, np.eye(n_in)); out = np.empty((L, n))
    for l in range(L):
        z = linear_step(st, W[l].astype(np.float64)); rec = {} if record is not None else None
        st = nonlin_step(z, o, rng, rec); out[l] = st.mu
        if record is not None: record[l] = rec
    return out
