"""Kikuchi algebras and adaptive final-mean states: the constructions of the 10 October 2026 note, implemented.

Every function is vectorised over final neurons (leading axis = neuron) unless stated otherwise.
Section numbers refer to the note.

  Sec 4   moment_interval              t_+ <= m <= (t+s)/2, Delta, the rank-two cut relaxation value
  Sec 5   restart                      anchor c, filtered state nu = (Z-c)^2 mu / v, bounded effect g_c, m = t_+ + Delta E_nu g_c
  Sec 6   quartic                      squared-shell defect d4, upper-gap bound d4/(4 s^3), effect s^2/(|Z|+s)^2
  Sec 7-8 spectral                     even (upper U_k) and odd (lower L_k) cyclic hierarchies from moments, residual R_k
  Sec 9   pair_overlap, three_lengths  m = E P - kappa E_nu 1/max(P,Q);  E relu(U_+ - V_+) closed form
  Sec 10  cell_interval, injection     conditional brackets from cell moments; sign-channel injection identity
  Sec 3   gate_tower                   exact up/down operators on a finite gate law (small N)
  Sec 11  kikuchi_laplacian, sparsify  all-level Kikuchi Laplacian, sampled sparsifier, common-congruence check
  Sec 12  schur                        coarse solve, residual identity, energy/trace certificates, enrichment
  Sec 13  observable_code_bound        RMS bound from constraint energy
"""
import numpy as np
from itertools import combinations

# ----------------------------------------------------------------------------------------------- Sec 4
def moment_interval(t, q):
    """First interval from t = E Z, q = E Z^2: [t_+, (t+s)/2], half-gap Delta = (s - |t|)/2."""
    t = np.asarray(t, float); s = np.sqrt(np.maximum(q, 0.0))
    return np.maximum(t, 0.0), 0.5 * (t + s), 0.5 * (s - np.abs(t))


def cut_relaxation_value(t, q):
    """sup over positive contractions 0 <= E <= I of Re <1, E Z> = top eigenvalue of K_Z = (|1><Z| + |Z><1|)/2."""
    t = np.asarray(t, float); s = np.sqrt(q)
    K = np.zeros(np.shape(t) + (2, 2))
    # orthonormal basis e1 = 1, e2 = (Z - t)/sqrt(q - t^2):  Z = t e1 + sqrt(q - t^2) e2
    b = np.sqrt(np.maximum(q - t * t, 0.0))
    K[..., 0, 0] = t; K[..., 0, 1] = K[..., 1, 0] = 0.5 * b
    return np.linalg.eigvalsh(K)[..., -1]


# ----------------------------------------------------------------------------------------------- Sec 5
def restart(t, q):
    """Anchor c = sign(t) s, normaliser v = E (Z-c)^2 = 4 s Delta, and the bounded effect g_c (a callable on z)."""
    t = np.asarray(t, float); s = np.sqrt(q); c = np.where(t >= 0, s, -s)
    lo, hi, D = moment_interval(t, q); v = 4.0 * s * D

    def g(z):
        z = np.asarray(z, float)
        R = np.where(c > 0, np.maximum(-z, 0.0), np.maximum(z, 0.0))
        den = (z - c) ** 2
        return np.where(den > 0, 4.0 * s * R / np.where(den > 0, den, 1.0), 0.0)
    return dict(c=c, v=v, Delta=D, lo=lo, g=g)


# ----------------------------------------------------------------------------------------------- Sec 6
def quartic(q, r):
    """d4 = r - q^2, upper-gap bound d4/(4 s^3) (eq. 11); effect s^2/(|Z|+s)^2 under nu4 = (Z^2-q)^2 mu / d4."""
    s = np.sqrt(q); d4 = r - q * q
    return dict(d4=d4, gap_bound=d4 / (4.0 * s ** 3), effect=lambda z: s ** 2 / (np.abs(z) + s) ** 2)


# ----------------------------------------------------------------------------------------------- Sec 7-8
def jacobi(mom, k, dps=60, rtol=1e-7):
    """Lanczos (Stieltjes) recurrence of the measure with ordinary moments mom[j] = int x^j, run in mpmath on the
    monomial coefficients of the orthonormal polynomials p_0, p_1, ...:
        x p_j = beta_{j-1} p_{j-1} + alpha_j p_j + beta_j p_{j+1}.
    Returns alpha[0..k-1] and beta[0..k-1] as mpf lists; beta[k-1] (the next residual norm) needs mom[0..2k], the
    alphas need mom[0..2k-1]. Missing moments shorten beta. A residual norm below rtol times the local scale is a
    breakdown (rtol reflects float64 moment inputs; exact inputs can use a smaller one)."""
    import mpmath as mp
    mp.mp.dps = dps; m = [mp.mpf(float(x)) if not isinstance(x, mp.mpf) else x for x in mom]
    ip = lambda u, v: mp.fsum(u[i] * v[j] * m[i + j] for i in range(len(u)) for j in range(len(v)))
    p_prev, p = [mp.mpf(0)], [1 / mp.sqrt(m[0])]; al, be = [], []
    for j in range(k):
        xp = [mp.mpf(0)] + p
        if 2 * j + 1 >= len(m): break
        a = ip(xp, p); al.append(a)
        r = [xp[i] - a * (p[i] if i < len(p) else 0) - (be[-1] * p_prev[i] if be and i < len(p_prev) else 0) for i in range(len(xp))]
        if 2 * j + 2 >= len(m): break
        b2 = ip(r, r); b = mp.sqrt(b2) if b2 > 0 else mp.mpf(0)
        if b <= rtol * (abs(a) + (be[-1] if be else 0)):     # breakdown (to the inputs' precision): exact termination
            break
        be.append(b)
        p_prev, p = p, [x / b for x in r]
    return al, be


def _eig_tridiag(al, be, k):
    """Eigen-decomposition of the k x k Jacobi matrix (mpmath, symmetric)."""
    import mpmath as mp
    J = mp.matrix(k, k)
    for i in range(k):
        J[i, i] = al[i]
        if i + 1 < k: J[i, i + 1] = J[i + 1, i] = be[i]
    E, Q = mp.eigsy(J)
    return [E[i] for i in range(k)], Q


def spectral(t, M, k, dps=60):
    """Even (upper U_k, Thm 7.1) and odd (lower L_k, Thm 7.2) cyclic hierarchies of one final coordinate from its
    ordinary moments M[p] = E Z^p (p = 0..4k; odd p are not used), with the residual certificate R_k (eqs. 15-16).
    Even measure: lambda = Z^2 under mu (moments M[2j]); odd measure: (lambda/q) times it (moments M[2j+2]/q).
    Returns floats U, L, R, the maintained interval [max(t_+, L, U - R), U], and the Gauss nodes and weights."""
    import mpmath as mp
    mp.mp.dps = dps; F = lambda x: x if isinstance(x, mp.mpf) else mp.mpf(float(x))     # mpf inputs keep full precision
    me = [F(M[2 * j]) for j in range(2 * k + 1) if 2 * j < len(M)]
    al, be = jacobi(me, k, dps); k = min(k, len(al))                      # a Lanczos breakdown means exact termination
    lam, V = _eig_tridiag(al, be, k)
    tiny = mp.mpf(10) ** (-dps // 2); sq = [mp.sqrt(max(x, tiny)) for x in lam]
    U = (F(t) + mp.fsum(V[0, a] ** 2 * sq[a] for a in range(k))) / 2
    R = mp.mpf(0) if len(be) < k and len(me) >= 2 * k + 1 else mp.mpf("nan")   # breakdown: B V_k in V_k, exact
    if len(be) >= k:
        h = [V[k - 1, a] * V[0, a] for a in range(k)]
        R = be[k - 1] ** 2 / 2 * mp.fsum(h[a] * h[b] / (sq[a] * sq[b] * (sq[a] + sq[b])) for a in range(k) for b in range(k))
    q = F(M[2]); mo = [F(M[2 * j + 2]) / q for j in range(2 * k) if 2 * j + 2 < len(M)]
    al2, be2 = jacobi(mo, k, dps); k2 = min(k, len(al2)); lam2, V2 = _eig_tridiag(al2, be2 + [mp.mpf(0)] * k2, k2)
    Lo = (F(t) + q * mp.fsum(V2[0, a] ** 2 / mp.sqrt(max(lam2[a], tiny)) for a in range(k2))) / 2
    U, Lo, R = float(U), float(Lo), float(R)
    return dict(U=U, L=Lo, R=R, lo=max(max(float(t), 0.0), Lo, U - R if np.isfinite(R) else -np.inf), hi=U,
                nodes=np.array([float(x) for x in lam]), weights=np.array([float(V[0, a] ** 2) for a in range(k)]))


def spectral_many(t, Mp, k):
    """t (n,), Mp (n, P) with Mp[:, p] = E Z^p. Returns arrays U, L, R, lo, hi (n,)."""
    out = [spectral(t[i], Mp[i], k) for i in range(len(t))]
    return {key: np.array([o[key] for o in out]) for key in ("U", "L", "R", "lo", "hi")}


# ----------------------------------------------------------------------------------------------- Sec 9
def three_lengths(u, v):
    """E ReLU((u.X)_+ - (v.X)_+) for X ~ N(0, I) (eq. 20)."""
    nu, nv, nd = np.linalg.norm(u), np.linalg.norm(v), np.linalg.norm(np.asarray(u) - np.asarray(v))
    return (nu - nv + nd) / (2 * np.sqrt(2 * np.pi))


def pair_split(w):
    """w = p - r with p, r >= 0 of disjoint support."""
    return np.maximum(w, 0.0), np.maximum(-w, 0.0)


def pair_overlap(prob, P, Q):
    """Prop. 9.1 on a finite (or empirical) law with atom probabilities prob and values P, Q >= 0:
    m = E P - kappa E_nu 1/max(P, Q), nu = P Q mu / kappa. Returns m, E P, kappa, E_nu 1/max, E_nu 1/max^2."""
    prob, P, Q = map(np.asarray, (prob, P, Q)); PQ = P * Q; kappa = np.sum(prob * PQ)
    mx = np.maximum(P, Q); on = PQ > 0; nu = np.where(on, prob * PQ, 0) / max(kappa, 1e-300)
    inv = np.where(on, 1 / np.where(on, mx, 1), 0)
    return dict(m=np.sum(prob * np.maximum(P - Q, 0)), EP=np.sum(prob * P), kappa=kappa,
                Enu_inv=np.sum(nu * inv), Enu_inv2=np.sum(nu * inv ** 2))


def pair_chain(prob, A, p, r):
    """The pair up/down chain of Sec. 9 on a finite law: atoms A (n_atoms, d) >= 0 with probabilities prob,
    final row split w = p - r. Pair states (a, b) with p_a r_b M_ab > 0; incidence Pi(a, b, A) = p_a r_b A_a A_b mu / kappa;
    S f(A) = E[f(a, b) | A]; K = S* S. Returns pi, K (on pair states), and the pair list."""
    prob, A = np.asarray(prob), np.asarray(A); P = A @ p; Q = A @ r; kappa = np.sum(prob * P * Q)
    pairs = [(a, b) for a in range(A.shape[1]) for b in range(A.shape[1]) if p[a] * r[b] * np.sum(prob * A[:, a] * A[:, b]) > 0]
    Pi = np.array([[prob[x] * p[a] * r[b] * A[x, a] * A[x, b] / kappa for x in range(len(prob))] for a, b in pairs])   # (pairs, atoms)
    pi = Pi.sum(1); nu = Pi.sum(0)                                    # pair marginal and state marginal (= P Q mu / kappa)
    cond = Pi / np.where(nu > 0, nu, 1)[None, :]                      # P(pair | atom)
    K = (Pi @ cond.T) / pi[:, None]                                   # K_{(ab),(a'b')} = sum_A Pi(ab | ...) ...
    return pi, K, pairs, nu, cond


def two_sign_heat_bath(c):
    """(E1 + E2)/2 on L^2 of two +-1 signs with correlation c (uniform marginals): eigenvalues 1, (1+c)/2, (1-c)/2, 0."""
    st = [(1, 1), (1, -1), (-1, 1), (-1, -1)]; pr = np.array([(1 + c) / 4, (1 - c) / 4, (1 - c) / 4, (1 + c) / 4])
    E = []
    for i in (0, 1):
        Ei = np.zeros((4, 4))
        for x, sx in enumerate(st):
            same = [y for y, sy in enumerate(st) if sy[i] == sx[i]]; w = pr[same] / pr[same].sum()
            Ei[x, same] = w
        E.append(Ei)
    T = (E[0] + E[1]) / 2; D = np.diag(np.sqrt(pr)); Ts = D @ T @ np.linalg.inv(D)      # symmetrise in L^2(pr)
    return np.sort(np.linalg.eigvals(Ts).real)[::-1]


# ----------------------------------------------------------------------------------------------- Sec 10
def cell_interval(t, p_cells, t_cells, q_cells):
    """Conditional brackets from unnormalised cell moments (cells on the last axis): lower sum (t_C)_+,
    upper (t + sum sqrt(p_C q_C))/2 (eq. 13), and the per-cell width contributions g_C = (sqrt(p_C q_C) - |t_C|)/2."""
    lo = np.sum(np.maximum(t_cells, 0.0), axis=-1)
    hi = 0.5 * (t + np.sum(np.sqrt(np.maximum(p_cells * q_cells, 0.0)), axis=-1))
    g = 0.5 * (np.sqrt(np.maximum(p_cells * q_cells, 0.0)) - np.abs(t_cells))
    return lo, hi, g


def conditional_brackets(prob, Z, labels):
    """Thm 10.1 on a finite law: for a field generated by a partition (atom labels), l = E (E[Z|F])_+ and
    u = (E Z + E sqrt(E[Z^2|F]))/2. Z may be (atoms,) or (atoms, n) for vector final coordinates."""
    prob = np.asarray(prob); Z = np.asarray(Z, float); lab = np.unique(labels, return_inverse=True)[1]
    pc = np.bincount(lab, prob); Zs = Z if Z.ndim == 2 else Z[:, None]
    t = np.stack([np.bincount(lab, prob * Zs[:, j]) for j in range(Zs.shape[1])], 1)
    q = np.stack([np.bincount(lab, prob * Zs[:, j] ** 2) for j in range(Zs.shape[1])], 1)
    lo = np.sum(np.maximum(t, 0), 0); hi = 0.5 * (np.sum(t, 0) + np.sum(np.sqrt(pc[:, None] * q), 0))
    return (lo, hi) if Z.ndim == 2 else (lo[0], hi[0])


def injection(z, h):
    """E_eps (z + eps h)_+ - z_+ = (|h| - |z|)_+ / 2 for a symmetric sign eps (eq. 23)."""
    return 0.5 * np.maximum(np.abs(h) - np.abs(z), 0.0)


def rms_bound(lo, hi):
    """Midpoint RMS certificate (eq. 14/21)."""
    return np.sqrt(np.sum((hi - lo) ** 2)) / (2 * np.sqrt(len(lo)))


# ----------------------------------------------------------------------------------------------- Sec 3
def gate_tower(prob, F, k):
    """Exact up/down operators on a finite gate law. prob: dict {bit-tuple of length N: probability}; F: dict
    {bit-tuple: final observable value (float)}. Returns f_S = E[F | G_S] for |S| = k, k+1, the matrices of U_k and
    D_{k+1} in the averaged L^2 inner products, and the checks of Theorem 3.1."""
    N = len(next(iter(prob)))
    configs = list(prob); P = np.array([prob[c] for c in configs]); Fv = np.array([F[c] for c in configs])

    def atoms(S):
        """atoms of F_S: map each config to its restriction b = c|_S; returns (labels, list of b)."""
        lab = {}; idx = []
        for c in configs:
            b = tuple(c[i] for i in S); idx.append(lab.setdefault(b, len(lab)))
        return np.array(idx), len(lab)

    def cond(S, vals):
        idx, na = atoms(S); num = np.bincount(idx, P * vals, na); den = np.bincount(idx, P, na)
        return (num / np.where(den > 0, den, 1))[idx]           # as a function on configs
    sets_k = list(combinations(range(N), k)); sets_k1 = list(combinations(range(N), k + 1))
    fk = {S: cond(S, Fv) for S in sets_k}; fk1 = {T: cond(T, Fv) for T in sets_k1}
    # up: (U_k f)_T = mean_{i in T} f_{T\{i}} ; down: (D h)_S = mean_{i notin S} E[h_{S+i} | F_S]
    Uf = {T: np.mean([fk[tuple(x for x in T if x != i)] for i in T], axis=0) for T in sets_k1}
    Dh = {S: np.mean([cond(S, fk1[tuple(sorted(S + (i,)))]) for i in range(N) if i not in S], axis=0) for S in sets_k}
    ip = lambda a, b, sets: np.mean([np.sum(P * a[S] * b[S]) for S in sets])
    # tower: D f_{k+1} = f_k ; adjointness <U f, h> = <f, D h> on random test vectors
    rng = np.random.default_rng(0)
    fr = {S: cond(S, rng.standard_normal(len(configs))) for S in sets_k}
    hr = {T: cond(T, rng.standard_normal(len(configs))) for T in sets_k1}
    Ufr = {T: np.mean([fr[tuple(x for x in T if x != i)] for i in T], axis=0) for T in sets_k1}
    Dhr = {S: np.mean([cond(S, hr[tuple(sorted(S + (i,)))]) for i in range(N) if i not in S], axis=0) for S in sets_k}
    adj = abs(ip(Ufr, hr, sets_k1) - ip(fr, Dhr, sets_k))
    tower = max(np.max(np.abs(Dh[S] - fk[S])) for S in sets_k)
    norm_id = (ip(fk1, fk1, sets_k1) - ip(fk, fk, sets_k)) - np.mean(
        [np.sum(P * (fk1[tuple(sorted(S + (i,)))] - fk[S]) ** 2) for S in sets_k for i in range(N) if i not in S])
    return dict(adjoint_err=adj, tower_err=tower, norm_identity_err=abs(norm_id), f=fk)


# ----------------------------------------------------------------------------------------------- Sec 11
def kikuchi_laplacian(N, wedges):
    """L_G = sum w_ij (I - V_ij) on l^2(2^[N]) (all levels at once). wedges: dict {(i,j): w}."""
    dim = 1 << N; Lm = np.zeros((dim, dim))
    for (i, j), w in wedges.items():
        for S in range(dim):
            bi, bj = (S >> i) & 1, (S >> j) & 1
            T = S ^ ((bi ^ bj) << i) ^ ((bi ^ bj) << j)              # swap slots i, j
            Lm[S, S] += w; Lm[S, T] -= w
    return Lm


def sparsify(wedges, q, rng):
    """Sample q edges with probability w_e / W and reweight by W/q (the sampled exchange part)."""
    keys = list(wedges); w = np.array([wedges[e] for e in keys]); W = w.sum()
    picks = rng.choice(len(keys), size=q, p=w / W); out = {}
    for p in picks:
        out[keys[p]] = out.get(keys[p], 0.0) + W / q
    return out


def relative_spectrum(Lm, Ls, tol=1e-10):
    """Extreme generalised eigenvalues of Ls against Lm on the range of Lm: (1-eps) Lm <= Ls <= (1+eps) Lm."""
    lam, V = np.linalg.eigh(Lm); keep = lam > tol * lam.max()
    Pm = V[:, keep] / np.sqrt(lam[keep]); M = Pm.T @ Ls @ Pm
    ev = np.linalg.eigvalsh(M); return ev.min(), ev.max()


# ----------------------------------------------------------------------------------------------- Sec 12
def schur(H, b, Cm, P, Q, gamma):
    """Coarse solve on retained coordinates P, exact correction (24), and the certificates of Theorem 12.1."""
    A = H[np.ix_(P, P)]; E = H[np.ix_(P, Q)]; Fm = H[np.ix_(Q, P)]; D = H[np.ix_(Q, Q)]
    Ai = np.linalg.inv(A); yP = Cm[:, P] @ Ai @ b[P]
    r = b[Q] - Fm @ Ai @ b[P]; S = D - Fm @ Ai @ E; K = Cm[:, Q] - Cm[:, P] @ Ai @ E
    y = Cm @ np.linalg.solve(H, b)
    FP = np.zeros_like(H); FP[np.ix_(P, P)] = Ai
    EP = b @ b / gamma - b @ FP @ b; UP = Cm @ Cm.T / gamma - Cm @ FP @ Cm.T
    n = Cm.shape[0]
    return dict(identity_err=np.abs(y - yP - K @ np.linalg.solve(S, r)).max(), err2=np.sum((y - yP) ** 2) / n,
                cert_spec=EP * np.linalg.eigvalsh(UP).max() / n, cert_trace=EP * np.trace(UP) / n, EP=EP, UP=UP)


def schur_enrich(H, b, Cm, P, J, gamma):
    """Thm 12.1 enrichment by a hidden block J: S_J, r_J, K_J and the updates of y_P, E_P, U_P, together with the
    trace-certificate decrease (delta_J Tr U_P + (E_P - delta_J) Tr Delta_J)/n."""
    A = H[np.ix_(P, P)]; EJ = H[np.ix_(P, J)]; Ai = np.linalg.inv(A)
    SJ = H[np.ix_(J, J)] - EJ.T @ Ai @ EJ; rJ = b[J] - EJ.T @ Ai @ b[P]; KJ = Cm[:, J] - Cm[:, P] @ Ai @ EJ
    base = schur(H, b, Cm, P, [i for i in range(len(b)) if i not in P], gamma)
    yP = Cm[:, P] @ Ai @ b[P]; SJi = np.linalg.inv(SJ)
    dJ = rJ @ SJi @ rJ; DJ = KJ @ SJi @ KJ.T; n = Cm.shape[0]
    return dict(y=yP + KJ @ SJi @ rJ, E=base["EP"] - dJ, U=base["UP"] - DJ, deltaJ=dJ, DeltaJ=DJ,
                trace_decrease=(dJ * np.trace(base["UP"]) + (base["EP"] - dJ) * np.trace(DJ)) / n)


# ----------------------------------------------------------------------------------------------- Sec 13
def observable_code_bound(Kop, Os, rho):
    """RMS((Tr rho O_i) - c_i) <= 2 b sqrt(q(1-q)) + d q (eq. 26), with P = ker K, gamma = smallest positive eigenvalue."""
    lam, V = np.linalg.eigh(Kop); ker = lam < 1e-12
    Pm = V[:, ker] @ V[:, ker].T; Qm = np.eye(len(lam)) - Pm; gamma = lam[~ker].min()
    c = np.array([np.trace(Pm @ O @ Pm) / max(np.trace(Pm), 1e-300) for O in Os])
    Bs = [O - ci * np.eye(len(lam)) for O, ci in zip(Os, c)]
    b2 = np.mean([np.linalg.norm(Pm @ B @ Qm, 2) ** 2 for B in Bs]); d2 = np.mean([np.linalg.norm(Qm @ B @ Qm, 2) ** 2 for B in Bs])
    eps = np.trace(rho @ Kop); q = np.trace(rho @ Qm)
    actual = np.sqrt(np.mean([(np.trace(rho @ O) - ci) ** 2 for O, ci in zip(Os, c)]))
    bound = 2 * np.sqrt(b2) * np.sqrt(max(q * (1 - q), 0)) + np.sqrt(d2) * q
    return dict(actual=actual, bound=bound, bound_eps=2 * np.sqrt(b2) * np.sqrt(eps / gamma) + np.sqrt(d2) * eps / gamma, q=q, eps_over_gamma=eps / gamma)
