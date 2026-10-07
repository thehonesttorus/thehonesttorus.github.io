# Shared machinery: Mehler sums with shifted Hermite coefficients, pair cumulants of h from raw moments, and the exact
# first variations of the raw pair moments of h = relu(y) under the third / fourth cumulant slices of y
# (Gaussian integration by parts: dE[F] = 1/6 sum kappa_ijk E[d_ijk F], 1/24 sum kappa_ijkl E[d_ijkl F]).
import numpy as np
from math import factorial
from scipy.special import ndtr
phi = lambda x: np.exp(-0.5 * x * x) / np.sqrt(2 * np.pi)
def mehler(P, Q, R, al, be, K):
    n = R.shape[0]; out = np.zeros((n, n)); Rj = np.ones((n, n))
    for j in range(K + 1):
        if j > 0: Rj = Rj * R
        out += np.outer(P[j + al], Q[j + be]) * (Rj / factorial(j))
    return out
def cum_pair(m, e2, e11, e21, e12, e22, e3, e4):
    C3 = e21 - e2[:, None] * m[None, :] - 2 * m[:, None] * e11 + 2 * (m * m)[:, None] * m[None, :]
    v = e2 - m * m; cov = e11 - np.outer(m, m)
    X2Y2 = (e22 - 2 * m[None, :] * e21 - 2 * m[:, None] * e12 + 4 * np.outer(m, m) * e11
            + (m * m)[None, :] * e2[:, None] + (m * m)[:, None] * e2[None, :] - 3 * np.outer(m * m, m * m))
    C4 = X2Y2 - np.outer(v, v) - 2 * cov * cov
    k4 = np.diag(C4).copy(); np.fill_diagonal(C4, 0.0)
    return C3, C4, k4
class Layer:
    """Gaussian reference of h = relu(y) for y with true (mu, S), and its first variations."""
    def __init__(self, mu, S, A, B, K=14):
        self.mu, self.S, self.A, self.B, self.K = mu, S, A, B, K
        self.sig = np.sqrt(np.diag(S)); self.R = S / np.outer(self.sig, self.sig); self.a = mu / self.sig
        a, sig, Pa, fa = self.a, self.sig, ndtr(self.a), phi(self.a); self.Pa, self.fa = Pa, fa
        self.mG = A[0]; self.e2G = B[0]
        self.e3G = sig**3 * ((a**3 + 3 * a) * Pa + (a * a + 2) * fa); self.e4G = sig**4 * ((a**4 + 6 * a * a + 3) * Pa + (a**3 + 5 * a) * fa)
        self.e11G = mehler(A, A, self.R, 0, 0, K); self.e21G = mehler(B, A, self.R, 0, 0, K); self.e12G = self.e21G.T.copy(); self.e22G = mehler(B, B, self.R, 0, 0, K)
        for M_, dg in ((self.e11G, self.e2G), (self.e21G, self.e3G), (self.e12G, self.e3G), (self.e22G, self.e4G)): np.fill_diagonal(M_, dg)
        self.C3G, self.C4G, self.k4G = cum_pair(self.mG, self.e2G, self.e11G, self.e21G, self.e12G, self.e22G, self.e3G, self.e4G)
    def _perturbed(self, dm, de2, d11, d21, d12, d22, de3, de4):
        for M_, dg in ((d11, de2), (d21, de3), (d12, de3), (d22, de4)): np.fill_diagonal(M_, dg)
        C3p, C4p, k4p = cum_pair(self.mG + dm, self.e2G + de2, self.e11G + d11, self.e21G + d21, self.e12G + d12, self.e22G + d22, self.e3G + de3, self.e4G + de4)
        return C3p - self.C3G, C4p - self.C4G, k4p - self.k4G
    def var3(self, k3, K21):
        """first variation of (C3, C4off, k4diag) of h under kappa3 slices of y: k3[a] = kappa(y_a,y_a,y_a), K21[a,c] = kappa(y_a,y_a,y_c)"""
        A, B, R, K, sig, a, Pa, fa = self.A, self.B, self.R, self.K, self.sig, self.a, self.Pa, self.fa
        k3s = k3 / sig**3; K21s = K21 / np.outer(sig**2, sig); K12s = K21s.T
        dm = k3s * A[3] / 6; de2 = k3s * B[3] / 6; de3 = k3s * sig**3 * Pa; de4 = 4 * k3s * sig**4 * (a * Pa + fa)
        def dE(P, Q):
            return (k3s[:, None] * mehler(P, Q, R, 3, 0, K) + 3 * K21s * mehler(P, Q, R, 2, 1, K)
                    + 3 * K12s * mehler(P, Q, R, 1, 2, K) + k3s[None, :] * mehler(P, Q, R, 0, 3, K)) / 6
        return self._perturbed(dm, de2, dE(A, A), dE(B, A), dE(A, B), dE(B, B), de3, de4)
    def var4(self, k4, K31, K22):
        """first variation under kappa4 slices: k4[a] = kappa(y_a^4), K31[a,c] = kappa(y_a,y_a,y_a,y_c), K22[a,c] = kappa(y_a,y_a,y_c,y_c)"""
        A, B, R, K, sig, a, Pa, fa = self.A, self.B, self.R, self.K, self.sig, self.a, self.Pa, self.fa
        k4s = k4 / sig**4; K31s = K31 / np.outer(sig**3, sig); K22s = K22 / np.outer(sig**2, sig**2); K13s = K31s.T
        dm = k4s * A[4] / 24; de2 = k4s * B[4] / 24; de3 = k4s * sig**3 * fa / 4; de4 = k4s * sig**4 * Pa
        def dE(P, Q):
            return (k4s[:, None] * mehler(P, Q, R, 4, 0, K) + 4 * K31s * mehler(P, Q, R, 3, 1, K) + 6 * K22s * mehler(P, Q, R, 2, 2, K)
                    + 4 * K13s * mehler(P, Q, R, 1, 3, K) + k4s[None, :] * mehler(P, Q, R, 0, 4, K)) / 24
        return self._perturbed(dm, de2, dE(A, A), dE(B, A), dE(A, B), dE(B, B), de3, de4)
def y_cumulants(D, l):
    """true cumulant slices of y_l from the MC raw moments"""
    mu = D["s1y"][l]; E2 = D["S2y"][l]; S = E2 - np.outer(mu, mu); e2 = np.diag(E2); sig2 = np.diag(S)
    k3 = D["s3y"][l] - 3 * mu * e2 + 2 * mu**3
    K21 = D["M21y"][l] - mu[None, :] * e2[:, None] - 2 * mu[:, None] * E2 + 2 * (mu * mu)[:, None] * mu[None, :]
    X4 = D["s4y"][l] - 4 * mu * D["s3y"][l] + 6 * mu * mu * e2 - 3 * mu**4
    X31 = (D["M31y"][l] - mu[None, :] * D["s3y"][l][:, None] - 3 * mu[:, None] * D["M21y"][l] + 3 * np.outer(mu, mu) * e2[:, None]
           + 3 * (mu * mu)[:, None] * E2 - 3 * np.outer(mu**3, mu))
    X22 = (D["M22y"][l] - 2 * mu[None, :] * D["M21y"][l] - 2 * mu[:, None] * D["M21y"][l].T + 4 * np.outer(mu, mu) * E2
           + (mu * mu)[None, :] * e2[:, None] + (mu * mu)[:, None] * e2[None, :] - 3 * np.outer(mu * mu, mu * mu))
    k4 = X4 - 3 * sig2**2; K31 = X31 - 3 * sig2[:, None] * S; K22 = X22 - np.outer(sig2, sig2) - 2 * S * S
    return dict(mu=mu, S=S, k3=k3, K21=K21, k4=k4, K31=K31, K22=K22)
def h_cumulants(D, l):
    m = D["s1h"][l]; e2 = np.diag(D["S2h"][l]); e11 = D["S2h"][l]; e21 = D["M21h"][l]
    return cum_pair(m, e2, e11, e21, e21.T.copy(), D["M22h"][l], D["s3h"][l], D["s4h"][l])
def cum_from_raw(s1, S2, M21, M31, M22, s3, s4):
    """cumulant slices (k3, K21, k4, K31, K22) from raw pair moments E y, E yy', E y^2 y', E y^3 y', E y^2 y'^2, E y^3, E y^4"""
    mu = s1; S = S2 - np.outer(mu, mu); e2 = np.diag(S2); sig2 = np.diag(S)
    k3 = s3 - 3 * mu * e2 + 2 * mu**3
    K21 = M21 - mu[None, :] * e2[:, None] - 2 * mu[:, None] * S2 + 2 * (mu * mu)[:, None] * mu[None, :]
    X4 = s4 - 4 * mu * s3 + 6 * mu * mu * e2 - 3 * mu**4
    X31 = (M31 - mu[None, :] * s3[:, None] - 3 * mu[:, None] * M21 + 3 * np.outer(mu, mu) * e2[:, None]
           + 3 * (mu * mu)[:, None] * S2 - 3 * np.outer(mu**3, mu))
    X22 = (M22 - 2 * mu[None, :] * M21 - 2 * mu[:, None] * M21.T + 4 * np.outer(mu, mu) * S2
           + (mu * mu)[None, :] * e2[:, None] + (mu * mu)[:, None] * e2[None, :] - 3 * np.outer(mu * mu, mu * mu))
    return dict(mu=mu, S=S, k3=k3, K21=K21, k4=X4 - 3 * sig2**2, K31=X31 - 3 * sig2[:, None] * S, K22=X22 - np.outer(sig2, sig2) - 2 * S * S)
def gauss_raw(m, C):
    """raw pair moments of a Gaussian vector with mean m, covariance C (orders <= 4, patterns (p,q) with p+q<=4)"""
    c = np.diag(C); n = len(m)
    mm = np.outer(m, m); E2 = C + mm
    # E[y_a^2 y_c] = m_c (c_a + m_a^2) + 2 m_a C_ac
    M21 = m[None, :] * (c + m * m)[:, None] + 2 * m[:, None] * C
    # E[y_a^3 y_c]: central (X+m_a)^3 (Y+m_c): E[X^3 Y] = 3 c_a C_ac ; E[X^2 Y] = 0; E[X Y]=C; E[X^2]=c
    M31 = (3 * c[:, None] * C + 3 * m[None, :] * (m * c)[:, None] * 0 + 3 * (m * m)[:, None] * C + 3 * m[:, None] * m[None, :] * c[:, None]
           + (m**3)[:, None] * m[None, :] + 3 * (m * c)[:, None] * m[None, :] * 0)
    # recompute M31 carefully: E[(X+a)^3 (Y+b)] with a=m_a, b=m_c:
    # = E[X^3 Y] + b E[X^3] + 3a E[X^2 Y] + 3ab E[X^2] + 3a^2 E[XY] + 3a^2 b E[X] + a^3 E[Y] + a^3 b
    # = 3 c_a C_ac + 0 + 0 + 3 a b c_a + 3 a^2 C_ac + 0 + 0 + a^3 b
    M31 = 3 * c[:, None] * C + 3 * (m * c)[:, None] * m[None, :] + 3 * (m * m)[:, None] * C + (m**3)[:, None] * m[None, :]
    # E[(X+a)^2 (Y+b)^2] = E[X^2Y^2] + 2b E[X^2 Y] + b^2 E[X^2] + 2a E[X Y^2] + 4ab E[XY] + 2ab^2 E[X] + a^2 E[Y^2] + 2a^2 b E[Y] + a^2 b^2
    #                    = c_a c_c + 2 C^2 + b^2 c_a + 4 a b C + a^2 c_c + a^2 b^2
    M22 = np.outer(c, c) + 2 * C * C + (m * m)[None, :] * c[:, None] + 4 * mm * C + (m * m)[:, None] * c[None, :] + np.outer(m * m, m * m)
    s3 = m**3 + 3 * m * c; s4 = m**4 + 6 * m * m * c + 3 * c * c
    return dict(s1=m.copy(), S2=E2, M21=M21, M31=M31, M22=M22, s3=s3, s4=s4)
def sm_slices(mu, S, gam, eps=1e-4):
    """cumulant slices of the scale mixture z = G y, y ~ N(mu, S), E G^2 = 1, Var G^2 = gam, to first order in gam
    (E G^r = 1 + c_r gam: c_1 = -1/8, c_2 = 0, c_3 = 3/8, c_4 = 1), by linearising the exact moment map at gam = 0"""
    raw = gauss_raw(mu, S); cr = {1: -1 / 8, 2: 0.0, 3: 3 / 8, 4: 1.0}
    def at(g):
        f = {r: 1 + cr[r] * g for r in cr}
        return cum_from_raw(raw["s1"] * f[1], raw["S2"] * f[2], raw["M21"] * f[3], raw["M31"] * f[4], raw["M22"] * f[4], raw["s3"] * f[3], raw["s4"] * f[4])
    p, m0 = at(eps), at(0.0)
    return {k: (p[k] - m0[k]) * (gam / eps) for k in ("k3", "K21", "k4", "K31", "K22")}
