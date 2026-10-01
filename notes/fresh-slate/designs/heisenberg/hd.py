"""Heisenberg-Duhamel (HD) estimator, numpy prototype (float64, full tensors, n <= 256).

Reference chain nu_l: Gaussian (m, C) of pre-activations, exact one-step moments.
Duhamel corrections: third cumulant tensors of the local non-Gaussianity ("sources"),
transported by the decoupled-gate linear response for up to A ages, paired with the
Heisenberg observable through Stein/Edgeworth injection into the next ReLU's mean and
covariance, then carried by the reference chain.

Full-tensor emulation: the n^3 tensors here give exactly the numbers the pulled-back
(source-layer) n x n evaluation of DESIGN.md section 2 gives; they are used only because
they are simpler to write at n <= 256.
"""
import numpy as np
from math import factorial
from scipy.special import ndtr, owens_t

SQ2PI = np.sqrt(2 * np.pi)


def phi(x):
    return np.exp(-0.5 * x * x) / SQ2PI


def hermite_e(k, x):
    """Probabilists' Hermite He_k(x)."""
    h0, h1 = np.ones_like(x), x
    if k == 0:
        return h0
    for j in range(1, k):
        h0, h1 = h1, x * h1 - j * h0
    return h1


def bvn_cdf(h, k, r):
    """P(X > -h, Y > -k)?  No: standard bivariate normal CDF Phi2(h, k; r) = P(X<h, Y<k), via Owen's T."""
    h = np.asarray(h, float); k = np.asarray(k, float); r = np.asarray(r, float)
    rr = np.sqrt(np.clip(1 - r * r, 1e-300, None))
    with np.errstate(divide="ignore", invalid="ignore"):
        ah = (k - r * h) / (h * rr)
        ak = (h - r * k) / (k * rr)
    out = 0.5 * ndtr(h) + 0.5 * ndtr(k)
    # handle h==0 or k==0 via limits: T(0,a)=atan(a)/(2pi); use tiny perturbation instead
    hz = np.abs(h) < 1e-12
    kz = np.abs(k) < 1e-12
    h2 = np.where(hz, 1e-12, h); k2 = np.where(kz, 1e-12, k)
    ah = (k2 - r * h2) / (h2 * rr)
    ak = (h2 - r * k2) / (k2 * rr)
    out = 0.5 * ndtr(h2) + 0.5 * ndtr(k2) - owens_t(h2, ah) - owens_t(k2, ak)
    out = out - 0.5 * ((h2 * k2 < 0))
    return out


class Gauss:
    """Gaussian quantities of z ~ N(m, C) needed by the ReLU layer."""

    def __init__(self, m, C, K=8):
        self.m, self.C = m, C
        self.s = np.sqrt(np.diag(C))
        self.t = m / self.s
        t, s = self.t, self.s
        self.Phi = ndtr(t)
        self.ph = phi(t)
        self.Ea = m * self.Phi + s * self.ph
        Ea2 = (m * m + s * s) * self.Phi + m * s * self.ph
        Ea3 = (m ** 3 + 3 * m * s * s) * self.Phi + s * (m * m + 2 * s * s) * self.ph
        self.var_a = Ea2 - self.Ea ** 2
        self.k3_a = Ea3 - 3 * self.Ea * Ea2 + 2 * self.Ea ** 3
        self.rho = C / np.outer(s, s)
        # E[delta^{(j)}(z)] = s^{-j-1} He_j(-t) phi(t)
        self.Ed = [s ** (-j - 1) * hermite_e(j, -t) * self.ph for j in range(K + 2)]
        # Hermite coefficients of centred a: alpha_1 = s Phi, alpha_k = s He_{k-2}(-t) phi(t)
        self.K = K
        self.alpha = np.zeros((K + 1, len(m)))
        self.alpha[1] = s * self.Phi
        for k in range(2, K + 1):
            self.alpha[k] = s * hermite_e(k - 2, -t) * self.ph
        # Hermite coefficients of centred a^2 (for coincident-index cumulants)
        self.g2 = np.zeros((K + 1, len(m)))
        Ea_z = self.Ea  # E[ReLU]
        self.g2[1] = 2 * s * Ea_z - 2 * self.Ea * self.alpha[1]
        self.g2[2] = 2 * s * s * self.Phi - 2 * self.Ea * self.alpha[2]
        for k in range(3, K + 1):
            self.g2[k] = 2 * s * s * hermite_e(k - 3, -t) * self.ph - 2 * self.Ea * self.alpha[k]
        # centred a^2 has mean Var(a); its centred version has the same k>=1 coefficients.

    def cov_a(self):
        """Exact Cov(ReLU(z_p), ReLU(z_q)) for the bivariate Gaussian."""
        m, s, r = self.m, self.s, self.rho
        n = len(m)
        r = np.clip(r, -1 + 1e-12, 1 - 1e-12)
        t1 = self.t[:, None] * np.ones((1, n))
        t2 = self.t[None, :] * np.ones((n, 1))
        rr = np.sqrt(1 - r * r)
        P2 = bvn_cdf(t1, t2, r)
        u1 = (t1 - r * t2) / rr
        u2 = (t2 - r * t1) / rr
        E = (t1 * t2 + r) * P2 + t1 * phi(t2) * ndtr(u1) + t2 * phi(t1) * ndtr(u2) + rr * phi(t1) * phi(u2)
        E = E * np.outer(s, s)
        Cv = E - np.outer(self.Ea, self.Ea)
        np.fill_diagonal(Cv, self.var_a)
        return Cv

    def cond(self):
        """Conditional law of z_q given z_p = 0: mean mu[p,q], sd sg[p,q]."""
        C, m = self.C, self.m
        d = np.diag(C)
        beta = C / d[:, None]  # beta[p,q] = C_pq / C_pp
        mu = m[None, :] - beta * m[:, None]
        var = d[None, :] - C * C / d[:, None]
        sg = np.sqrt(np.clip(var, 1e-300, None))
        return beta, mu, sg

    def source_k3(self, Kt=4, diagrams="all"):
        """Third cumulant tensor of a = ReLU(z), z Gaussian (centred moments)."""
        n = len(self.m)
        al, rho = self.alpha, self.rho
        T = np.zeros((n, n, n))
        if diagrams in ("all", "star", "startri"):
            rp = [np.ones_like(rho)] + [rho ** i for i in range(1, Kt + 1)]
            for i in range(Kt + 1):
                for j in range(Kt + 1):
                    for k in range(Kt + 1):
                        a, b, c = i + j, i + k, j + k
                        if a < 1 or b < 1 or c < 1 or max(a, b, c) > self.K:
                            continue
                        if diagrams == "star" and not ((i, j, k) in ((1, 1, 0), (1, 0, 1), (0, 1, 1))):
                            continue
                        if diagrams == "startri" and not ((i, j, k) in ((1, 1, 0), (1, 0, 1), (0, 1, 1), (1, 1, 1))):
                            continue
                        coef = 1.0 / (factorial(i) * factorial(j) * factorial(k))
                        A = (al[a][:, None] * rp[i])  # p,q
                        Bm = (al[c][None, :] * rp[j])  # p,r  (alpha_r)
                        T += coef * (A[:, :, None] * Bm[:, None, :]) * (al[b][:, None] * rp[k])[None, :, :]
        # coincident indices: p=q!=r etc. via a^2 coefficients; p=q=r exact
        Kc = self.K
        rpw = [rho ** k for k in range(Kc + 1)]
        P = np.zeros((n, n))  # P[p, r] = E[ãp^2 ãr]
        for k in range(1, Kc + 1):
            P += self.g2[k][:, None] * al[k][None, :] * rpw[k] / factorial(k)
        idx = np.arange(n)
        T[idx, idx, :] = P
        T[idx, :, idx] = P
        T[:, idx, idx] = P.T
        T[idx, idx, idx] = self.k3_a
        return T


def tmode(T, W):
    """Contract all three legs of a symmetric tensor with W (x @ W convention): T'[a,b,c]=sum W_pa W_qb W_rc T_pqr."""
    n = T.shape[0]
    X = np.tensordot(T, W, axes=([2], [0]))          # p q c
    X = np.tensordot(X, W, axes=([1], [0]))          # p c b
    X = np.tensordot(X, W, axes=([0], [0]))          # c b a
    return X.transpose(2, 1, 0)


GH_X, GH_W = np.polynomial.hermite_e.hermegauss(200)
GH_W = GH_W / GH_W.sum()


def k4_single(G):
    """Exact fourth cumulant of ReLU(z_p), z_p ~ N(m_p, s_p^2), by quadrature (kink handled by splitting)."""
    out = np.empty(len(G.m))
    # E[X+^k] closed forms via recursion on truncated normal moments
    t, s, m = G.t, G.s, G.m
    # moments of X+ : use M_k = E[(m + sZ)^k 1{Z > -t}] via binomial + truncated std normal moments
    Phi, ph = G.Phi, G.ph
    # truncated moments T_j = E[Z^j 1{Z > -t}]: T0 = Phi, T1 = phi(t), T_j = (-t)^{j-1} phi(t) + (j-1) T_{j-2}
    T = [Phi, ph]
    for j in range(2, 5):
        T.append((-t) ** (j - 1) * ph + (j - 1) * T[j - 2])
    from math import comb
    M = [sum(comb(k, j) * m ** (k - j) * s ** j * T[j] for j in range(k + 1)) for k in range(5)]
    mu = M[1]
    c2 = M[2] - mu ** 2
    c4 = M[4] - 4 * mu * M[3] + 6 * mu ** 2 * M[2] - 3 * mu ** 4
    return c4 - 3 * c2 ** 2


def gen_k4_diag(G, W, cycle=True):
    """Diagonal fourth cumulant of z' = ReLU(z) W for Gaussian z: second-chaos diagrams (chain, cycle),
    the degree-3 star, and the exact single-site term."""
    al, rho = G.alpha, G.rho
    B = rho @ (al[1][:, None] * W)                 # (Sigma b)_p for every output i
    Qd = 0.5 * al[2][:, None] * W                  # Q_pp per output
    U = rho @ (Qd * B)
    chain = 48 * np.sum(Qd * B * U, axis=0)
    star = 4 * np.sum(al[3][:, None] * W * B ** 3, axis=0)
    cyc = np.zeros(W.shape[1])
    if cycle:
        for i in range(W.shape[1]):
            M = rho * Qd[:, i][None, :]
            M2 = M @ M
            cyc[i] = 48 * np.sum(M2 * M2.T)
    # single-site: replace the diagonal parts of the diagrams by the exact single-site cumulant
    a1, a2, a3 = al[1], al[2], al[3]
    quad_single = 3 * a2 ** 4 + 12 * a1 ** 2 * a2 ** 2 + 4 * a3 * a1 ** 3
    single = (k4_single(G) - quad_single)[:, None] * W ** 4
    return chain + star + cyc + single.sum(0)


def inject2(G, D, K4):
    """Second-order Edgeworth on the marginals: kappa_4 (K4[p]) and kappa_3^2 terms on mean and variance."""
    dEa = K4 * G.Ed[2] / 24.0 + D * D * G.Ed[4] / 72.0
    dEa2 = K4 * 2 * G.Ed[1] / 24.0 + D * D * 2 * G.Ed[3] / 72.0
    return dEa, dEa2


def inject(G, D, S):
    """Stein/Edgeworth kappa_3 injection: returns dEa (n,), dCov (n,n).
    D[p] = k3(z_p,z_p,z_p), S[p,q] = k3(z_p,z_p,z_q)."""
    n = len(G.m)
    dEa = D * G.Ed[1] / 6.0
    beta, mu, sg = G.cond()
    p0 = G.Ed[0]                                         # density of z_p at 0
    PhiC = ndtr(mu / sg)                                  # P(z_q>0 | z_p=0)
    EreluC = mu * PhiC + sg * phi(mu / sg)                # E[ReLU z_q | z_p=0]
    EdH = p0[:, None] * PhiC                              # E[delta(z_p) H(z_q)]
    pprime0 = (G.m / G.s ** 2) * p0                       # p'(0)
    Edpa = -(pprime0[:, None] * EreluC + p0[:, None] * beta * PhiC)   # E[delta'(z_p) ReLU(z_q)]
    dEaa = (D[:, None] * Edpa + 3 * S * EdH + 3 * S.T * EdH.T + D[None, :] * Edpa.T) / 6.0
    dC = dEaa - dEa[:, None] * G.Ea[None, :] - G.Ea[:, None] * dEa[None, :]
    dvar = D * p0 / 3.0 - 2 * G.Ea * dEa
    np.fill_diagonal(dC, dvar)
    return dEa, dC


def closure(Ws):
    """Reference Gaussian closure. Ws: (L, n, n) float64. Returns (L, n) means of a_l."""
    L, n, _ = Ws.shape
    m = np.zeros(n); C = Ws[0].T @ Ws[0]
    out = []
    for l in range(L):
        G = Gauss(m, C, K=2)
        out.append(G.Ea)
        if l + 1 < L:
            Ca = G.cov_a()
            m = G.Ea @ Ws[l + 1]
            C = Ws[l + 1].T @ Ca @ Ws[l + 1]
    return np.array(out)


def hd(Ws, A=None, diagrams="all", Kt=4, K=8, only_src=None, second=False, cycle=False):
    """HD-A estimator. A=None: all ages. Returns (L, n) means."""
    L, n, _ = Ws.shape
    m = np.zeros(n); C = Ws[0].T @ Ws[0]
    K4 = None
    K4s = None
    ages = []   # list of kappa3 tensors of z_l by age (0 = generated by the last arrow)
    out = []
    for l in range(L):
        G = Gauss(m, C, K=max(K, 21) if second == "full" else K)
        if ages:
            Kz = sum(ages)
            D = np.einsum("ppp->p", Kz).copy()
            S = np.einsum("ppq->pq", Kz).copy()
            if second == "full" and K4s is not None:
                dEa, dC = inject_full2(G, D, S, K4s[2], K4s[1], K4s[0])
            else:
                dEa, dC = inject(G, D, S)
            if second is True:
                d2, d2a = inject2(G, D, K4 if K4 is not None else np.zeros(n))
                dEa = dEa + d2
                dC[np.diag_indices(n)] += d2a - 2 * G.Ea * d2
        else:
            dEa, dC = np.zeros(n), np.zeros((n, n))
        Ea = G.Ea + dEa
        out.append(Ea)
        if l + 1 < L:
            W = Ws[l + 1]
            Ca = G.cov_a() + dC
            if only_src is None or only_src == l:
                src = G.source_k3(Kt=Kt, diagrams=diagrams)
                new = [tmode(src, W)]
            else:
                new = [np.zeros((n, n, n))]
            if A is None or A > 0:
                g = G.Phi
                for a, T in enumerate(ages):
                    if A is not None and a + 1 > A:
                        break
                    Tg = T * g[:, None, None] * g[None, :, None] * g[None, None, :]
                    new.append(tmode(Tg, W))
                if A is None and len(new) > 1:
                    new = [new[0], sum(new[1:])]  # keep memory bounded: age 0 and "older"
            ages = new if A is not None else new
            if second == "full":
                K4s = gen_k4_slices(G, W, cycle=cycle, pairs_exact=True)
            elif second and (only_src is None or only_src == l):
                K4 = gen_k4_diag(G, W)
            else:
                K4 = None
            m = Ea @ W
            C = W.T @ Ca @ W
    return np.array(out)


def ederiv(G, r):
    """E[a^{(r)}(z_p)] in z-units: r=0 ReLU mean, r=1 Phi, r>=2 E delta^{(r-2)}."""
    if r == 0:
        return G.Ea
    if r == 1:
        return G.Phi
    return G.Ed[r - 2]


def biv(G, i, j, K=14):
    """E[a^{(i)}(z_p) a^{(j)}(z_q)] for all pairs via the Mehler series (z-units).
    Needs G built with K >= i+K-2 etc.; diagonal is not meaningful (overwritten by callers)."""
    rho = G.rho; s = G.s
    out = np.zeros_like(rho)
    rk = np.ones_like(rho); fk = 1.0
    for k in range(K + 1):
        if k > 0:
            rk = rk * rho; fk *= k
        A = ederiv(G, i + k) * s ** k
        B = ederiv(G, j + k) * s ** k
        out += rk / fk * np.outer(A, B)
    return out


def inject_full2(G, D, S, K4d, K22, K31, K=14, k3sq=True):
    """Second-order Edgeworth on mean and the FULL covariance of a = ReLU(z).
    D: k3 diag, S[p,q] = k3(ppq), K4d: k4 diag, K22[p,q] = k4(ppqq), K31[p,q] = k4(pppq)."""
    from math import comb
    n = len(G.m)
    dEa = D * G.Ed[2 - 1] / 6 + K4d * G.Ed[2] / 24 + D * D * G.Ed[4] / 72
    E = lambda i, j: biv(G, i, j, K)
    c3 = {3: D[:, None] * np.ones((1, n)), 2: S, 1: S.T, 0: np.ones((n, 1)) * D[None, :]}
    # third order
    d = (c3[3] * E(3, 0) + 3 * c3[2] * E(2, 1) + 3 * c3[1] * E(1, 2) + c3[0] * E(0, 3)) / 6
    # fourth order
    c4 = {4: K4d[:, None] * np.ones((1, n)), 3: K31, 2: K22, 1: K31.T, 0: np.ones((n, 1)) * K4d[None, :]}
    d += sum(comb(4, k) * c4[k] * E(k, 4 - k) for k in range(5)) / 24
    # kappa3^2
    for k in (range(7) if k3sq else []):
        coef = 0
        for i in range(4):
            j = k - i
            if 0 <= j <= 3:
                coef = coef + comb(3, i) * comb(3, j) * c3[i] * c3[j]
        if not np.isscalar(coef):
            d += comb(6, k) / comb(6, k) * coef * E(k, 6 - k) / 72
    dC = d - dEa[:, None] * G.Ea[None, :] - G.Ea[:, None] * dEa[None, :] - np.outer(dEa, dEa)
    dvar = (D * 2 * G.Ed[0] / 6 + K4d * 2 * G.Ed[1] / 24 + D * D * 2 * G.Ed[3] / 72) - 2 * G.Ea * dEa - dEa ** 2
    np.fill_diagonal(dC, dvar)
    return dEa, dC


def gen_k4_slices(G, W, cycle=True, pairs_exact=True):
    """Generated joint fourth-cumulant slices of y = ReLU(z) W, z ~ N(m, C):
    K31[p,q] = k4(y_p,y_p,y_p,y_q), K22[p,q] = k4(y_p,y_p,y_q,y_q), K4d[p] = k4(y_p).
    Diagrams of the second-chaos CGF (chain, cycle) + degree-3 star + exact single-site."""
    al, Sg = G.alpha, G.rho
    n = W.shape[1]
    B = Sg @ (al[1][:, None] * W)
    Qd = 0.5 * al[2][:, None] * W
    V = Sg @ (Qd * B)
    a3 = al[3][:, None]
    T1 = (Qd * V).T @ B
    T2 = (B * V).T @ Qd
    K31 = 24 * (T1 + T2)
    K31 += (a3 * B ** 3).T @ W + 3 * (a3 * W * B ** 2).T @ B
    a1, a2 = al[1], al[2]
    quad_single = 3 * a2 ** 4 + 12 * a1 ** 2 * a2 ** 2 + 4 * al[3] * a1 ** 3
    cs = (k4_single(G) - quad_single)[:, None]
    K31 += (cs * W ** 3).T @ W
    # K22
    X1 = V.T @ (Qd * B)                      # bp S Qp S Qq S bq = sum_r V_rp Qd_rq B_rq
    K22 = 2 * ((W * B).T @ (a3 * B * B)) + 2 * ((a3 * B * B).T @ (W * B))   # star: 2 sum a3 (W_p B_p B_q^2 + W_q B_p^2 B_q)
    K22 += (cs * W ** 2).T @ W ** 2
    chain = 2 * X1 + 2 * X1.T                # 2 bpSQpSQqSbq + 2 bqSQqSQpSbp (transposes)
    if pairs_exact:
        P1 = np.empty((n, n)); P2 = np.empty((n, n))
        for p in range(n):
            BQ = B[:, p:p + 1] * Qd             # (B_p o Qd_q) columns q
            P1[p] = np.sum(BQ * (Sg @ BQ), axis=0)               # bp S Qq S Qq S bp
            QB = Qd[:, p:p + 1] * B             # (Qd_p o B_q)
            P2[p] = np.sum(BQ * (Sg @ QB), axis=0)               # bp S Qq S Qp S bq
        chain += P1 + P1.T + 2 * P2
    K22 += 8 * chain
    if cycle:
        for p in range(n):
            M = Sg * Qd[:, p][None, :]                # S Qp
            M2 = M @ M
            M3S = M2 @ M @ Sg                         # (S Qp)^3 S
            K31[p] += 48 * (np.diag(M3S)[:, None] * Qd).sum(0)
            # s^2 t^2: 4 tr(SQp SQp SQq SQq) + 2 tr(SQp SQq SQp SQq), /6, x48
            A = M2 @ Sg                               # S Qp S Qp S
            # tr(S Qp S Qp S Qq S Qq) = sum_{r,s} A_rs Qq_s S_sr Qq_r -> Qq^T (A o S) Qq
            t1 = np.sum(Qd * ((A * Sg.T) @ Qd), axis=0)
            Bm = M @ Sg                               # S Qp S
            t2 = np.sum(Qd * ((Bm * Bm.T) @ Qd), axis=0)
            K22[p] += 48 * (4 * t1 + 2 * t2) / 6
    return K31, K22, np.diag(K31).copy()


def tmode4(T, W):
    for _ in range(4):
        T = np.tensordot(T, W, axes=([0], [0]))   # rotates legs; after 4 passes order restored
    return T


def hd2(Ws, A4=None, k4src="single", Kt=4):
    """HD with all-age kappa_3 (full series) AND all-age kappa_4 tensor transported with decoupled gates;
    kappa_4 source = exact single-site kappa_4(a_r) (diagonal tensor). Full second-order injection.
    Full n^4 tensors: prototype only (n <= 64)."""
    L, n, _ = Ws.shape
    m = np.zeros(n); C = Ws[0].T @ Ws[0]
    K3 = None; K4 = None; out = []
    for l in range(L):
        G = Gauss(m, C, K=21)
        if K3 is not None:
            D = np.einsum("ppp->p", K3).copy(); S = np.einsum("ppq->pq", K3).copy()
            K4d = np.einsum("pppp->p", K4).copy(); K31 = np.einsum("pppq->pq", K4).copy(); K22 = np.einsum("ppqq->pq", K4).copy()
            dEa, dC = inject_full2(G, D, S, K4d, K22, K31)
        else:
            dEa, dC = np.zeros(n), np.zeros((n, n))
        Ea = G.Ea + dEa; out.append(Ea)
        if l + 1 < L:
            W = Ws[l + 1]; g = G.Phi
            src3 = G.source_k3(Kt=Kt)
            new3 = src3 if K3 is None else src3 + K3 * g[:, None, None] * g[None, :, None] * g[None, None, :]
            K3 = tmode(new3, W)
            c = k4_single(G)
            src4 = np.einsum("r,ra,rb,rc,rd->abcd", c, W, W, W, W, optimize=True)
            if K4 is None:
                K4 = src4
            else:
                g4 = K4 * g[:, None, None, None] * g[None, :, None, None] * g[None, None, :, None] * g[None, None, None, :]
                K4 = src4 + tmode4(g4, W)
            Ca = G.cov_a() + dC
            m = Ea @ W; C = W.T @ Ca @ W
    return np.array(out)
