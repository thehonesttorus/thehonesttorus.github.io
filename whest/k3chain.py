"""Third-cumulant chain in CP (leg) form, derived from the Wiener-chaos / Hermite picture.

Post-activation third cumulant born at layer b (h = relu(z), z ~ N(mu, C) to first order in the off-diagonal
correlations, exact in the gate margins alpha = mu/sigma):
    T_b = sum_m w2_m Sym3(e_m x a_m x a_m)  +  sum_m 3 s_m Sym3(e_m x e_m x a_m)  +  sum_m d_m e_m x e_m x e_m,
with w2 = phi(alpha)/sigma (the kink strength E relu''), the arm a_m = D(Phi) C e_m (= Cov(h, z_m) to first order),
s_m the exact first-order (2,1) slice coefficient and d_m the exact 1-D third cumulant correction.
Transport: through W the legs multiply, through relu the legs are gated by Phi (product gate). Readouts at layer l:
    D3_i   = kappa3(z_i)        = sum_m w2 P A^2 + s P^2 A + d P^3   (per source)
    D21_ij = kappa3(z_i,z_i,z_j) = [ (2/3) P.A w2 + P.P s ] A^T  +  [ (1/3) A.A w2 + 2 P.A s + P.P d ] P^T.
Mean readout (Edgeworth): E relu(z) = sigma d0(alpha) - kappa3 alpha phi/(6 sigma^2) + kappa4 (alpha^2-1) phi/(24 sigma^3).
Second moment: E relu^2 += kappa3 phi/(3 sigma); covariance: Cov(h_i,h_j) += (D21_ij phi_i Phi_j/sigma_i + sym)/2.
Options (dict): fold (0/1: newborn arm += D(phi/sigma) D21 / 2), hub (0/1: D21 readouts), k4 ('none' | 'path':
the second-chaos path class 12 diag(Y C^-1 Y^T) with Y the w2-weighted hub product, plus the third-chaos star 4 sum
c3 P A^3), K (Hermite order), window (max source age, 0 = all).
"""
import numpy as np
from scipy.special import ndtr
from .relu_gauss import phi, hermite_relu


def k3_exact(alpha):
    """third cumulant of relu(u + alpha), u ~ N(0,1)."""
    P = ndtr(alpha); p = phi(alpha)
    m1 = p + alpha * P; m2 = (1 + alpha ** 2) * P + alpha * p; m3 = (alpha ** 2 + 2) * p + (3 * alpha + alpha ** 3) * P
    return m3 - 3 * m2 * m1 + 2 * m1 ** 3


class Source:
    __slots__ = ("P", "A", "w2", "s", "d", "born", "c3")
    def __init__(self, P, A, w2, s, d, born, c3):
        self.P, self.A, self.w2, self.s, self.d, self.born, self.c3 = P, A, w2, s, d, born, c3


def k3_chain(W, opts=None, record=None):
    o = dict(fold=0, hub=1, k4="none", K=8, window=0, eps=0.01)
    if opts: o.update(opts)
    L, n, n_in = W.shape
    m = np.zeros(n_in); Kh = np.eye(n_in); sources = []; out = np.empty((L, n)); flops = 0.0
    for l in range(L):
        Wl = W[l].astype(np.float64); last = (l == L - 1)
        mu = Wl @ m; C = Wl @ Kh @ Wl.T; flops += 4 * n ** 3
        for s in sources:
            s.P = Wl @ s.P; s.A = Wl @ s.A; flops += 4 * n ** 3
        var = np.clip(np.diag(C), 1e-30, None); sigma = np.sqrt(var); alpha = mu / sigma
        Phi = ndtr(alpha); ph = phi(alpha)
        D3 = np.zeros(n); D21 = np.zeros((n, n)); Y = np.zeros((n, n)); star4 = np.zeros(n)
        for s in sources:
            PA = s.P * s.A
            D3 += (PA * s.A) @ s.w2 + (s.P * s.P * s.A) @ s.s + (s.P ** 3) @ s.d
            if o["k4"] == "path":
                Yc = (PA * s.w2[None, :]) @ s.A.T; Y += Yc; flops += 2 * n ** 3
                star4 += (s.P * s.A ** 3) @ s.c3
            if o["hub"] and not last:
                LA = (2.0 / 3.0) * PA * s.w2[None, :] + (s.P * s.P) * s.s[None, :]
                LP = (1.0 / 3.0) * (s.A * s.A) * s.w2[None, :] + 2.0 * PA * s.s[None, :] + (s.P * s.P) * s.d[None, :]
                D21 += LA @ s.A.T + LP @ s.P.T; flops += 4 * n ** 3
        d = hermite_relu(alpha, o["K"])
        K4 = np.zeros(n)
        if o["k4"] == "path" and sources:
            Creg = C + o["eps"] * np.mean(var) * np.eye(n)
            Z = np.linalg.solve(Creg, Y.T); flops += (2.0 / 3.0) * n ** 3 + 2 * n ** 3
            K4 = 12.0 * np.einsum("ij,ji->i", Y, Z) + 4.0 * star4
            K4 -= np.mean(K4)
        m_new = sigma * d[0] - D3 * alpha * ph / (6 * var) + K4 * (alpha ** 2 - 1) * ph / (24 * sigma ** 3)
        out[l] = m_new
        if record is not None:
            record[l] = dict(mu=mu, var=var, alpha=alpha, D3=D3.copy(), D21=D21.copy(), m=m_new.copy(), K4=K4.copy())
        if last:
            break
        second = var * ((1 + alpha ** 2) * Phi + alpha * ph) + D3 * ph / (3 * sigma)
        rho = C / np.outer(sigma, sigma); np.fill_diagonal(rho, 0.0)
        Kh_new = np.zeros_like(C); rk = np.ones_like(rho); fact = 1.0
        for k in range(1, o["K"] + 1):
            rk = rk * rho; fact *= k; Kh_new += np.outer(d[k], d[k]) * rk / fact
        Kh_new *= np.outer(sigma, sigma)
        if o["hub"]:
            g = ph / sigma
            Kh_new += 0.5 * (D21 * g[:, None] * Phi[None, :] + D21.T * Phi[:, None] * g[None, :])
        np.fill_diagonal(Kh_new, second - m_new ** 2)
        # gate the carried sources into post-activation space
        for s in sources:
            s.P *= Phi[:, None]; s.A *= Phi[:, None]
        if o["window"]:
            sources = [s for s in sources if l + 1 - s.born <= o["window"]]
        # the newborn source
        w2 = ph / sigma; mg = sigma * d[0]
        A = Phi[:, None] * C
        if o["fold"] and o["hub"]:
            A = A + 0.5 * (ph / sigma)[:, None] * D21
        e = 2 * mg * (1 - Phi)
        s_vec = e - (2.0 / 3.0) * w2 * Phi ** 2 * var
        d_vec = sigma ** 3 * k3_exact(alpha) - w2 * (Phi * var) ** 2 - s_vec * Phi * var
        c3 = -alpha * ph / var          # E relu''' = -alpha phi / sigma^2
        sources.append(Source(np.eye(n), A, w2, s_vec, d_vec, l, c3))
        m, Kh = m_new, Kh_new
    return out, flops
