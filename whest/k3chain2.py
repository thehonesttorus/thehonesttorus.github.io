"""Third-cumulant chain, version 2: CP legs P (centre), A (covariance arm), B (third-cumulant arm: the GC terms),
R (slice/diagonal residual, exact at birth: all Hermite orders of the (2,1) slice, the exact diagonal, and the exact
ReLU gates of the transported sources' slices and diagonals). Terms (m = birth neuron, legs = columns):
    star:   sum_m w2_m  (e x a x a + a x e x a + a x a x e)_m            [w2 = E relu'' = phi/sigma]
    GC:     sum_m       sum_{S_3} (e x a x b)_m                          [b_m = w3_m Phi.D21[m,:] + w2_m/2 w2.D21[:,m]]
    res:    sum_m       (e x e x r + e x r x e + r x e x e)_m
Readouts in z-space:  D3 = 3 w2 P A^2 + 6 P A B + 3 P^2 R  (summed over m);
    D21_ij = sum_m w2 (2 P_i A_i A_j + A_i^2 P_j) + 2 (P_i A_i B_j + P_i B_i A_j + A_i B_i P_j) + (P_i^2 R_j + 2 P_i R_i P_j).
Dense legs throughout (an accuracy instrument; cost = 8 n^3 per source-layer).
"""
import numpy as np
from scipy.special import ndtr
from .relu_gauss import phi, hermite_relu
from .k3chain import k3_exact


def hermite_neg(x, J):
    """He_j(-x) for j = 0..J, shape (J+1, *x.shape)."""
    H = np.empty((J + 1,) + x.shape); H[0] = 1.0
    if J >= 1: H[1] = -x
    for j in range(2, J + 1):
        H[j] = -x * H[j - 1] - (j - 1) * H[j - 2]
    return H


def relu2_hermite(alpha, sigma, K):
    """Hermite coefficients (w.r.t. z ~ N(mu, sigma^2)) of relu(z)^2:  k = 0..K.
    F_0 = sigma^2 m2(alpha); F_1 = 2 sigma d0; F_2 = 2 Phi; F_k = 2 sigma^{2-k} He_{k-3}(-alpha) phi(alpha), k >= 3."""
    P = ndtr(alpha); p = phi(alpha); d0 = p + alpha * P
    F = np.empty((K + 1,) + alpha.shape)
    F[0] = sigma ** 2 * ((1 + alpha ** 2) * P + alpha * p)
    if K >= 1: F[1] = 2 * sigma * d0
    if K >= 2: F[2] = 2 * P
    if K >= 3:
        H = hermite_neg(alpha, K - 3)
        for k in range(3, K + 1):
            F[k] = 2 * sigma ** (2 - k) * H[k - 3] * p
    return F


def relu_hermite_z(alpha, sigma, K):
    """Hermite coefficients w.r.t. z ~ N(mu, sigma^2) of relu(z): d_k(alpha) sigma^{1-k}."""
    d = hermite_relu(alpha, K)
    return d * sigma[None, :] ** (1 - np.arange(K + 1))[:, None]


class Src:
    __slots__ = ("P", "A", "B", "R", "w2", "born")
    def __init__(self, P, A, B, R, w2, born):
        self.P, self.A, self.B, self.R, self.w2, self.born = P, A, B, R, w2, born


def k3_chain2(W, opts=None, record=None):
    o = dict(K=8, KS=8, gc=1, res=1, gate_res=1, k4="none", eps=0.01, window=0, hub=1, k4diag=None, k4var=1, k22born=0, k31born=0, pert_off=0.0, pert_var=0.0, d3scale=1.0)
    if opts: o.update(opts)
    L, n, n_in = W.shape
    m = np.zeros(n_in); Kh = np.eye(n_in); srcs = []; out = np.empty((L, n)); flops = 0.0
    K22z = None; K31z = None
    for l in range(L):
        Wl = W[l].astype(np.float64); last = (l == L - 1)
        mu = Wl @ m; C = Wl @ Kh @ Wl.T; flops += 4 * n ** 3
        for s in srcs:
            s.P = Wl @ s.P; s.A = Wl @ s.A; s.B = Wl @ s.B; s.R = Wl @ s.R; flops += 8 * n ** 3
        var = np.clip(np.diag(C), 1e-30, None); sigma = np.sqrt(var); alpha = mu / sigma
        Phi = ndtr(alpha); ph = phi(alpha)
        D3 = np.zeros(n); D21 = np.zeros((n, n)); Y = np.zeros((n, n)); star4 = np.zeros(n)
        for s in srcs:
            PA = s.P * s.A
            D3 += 3.0 * ((PA * s.A) @ s.w2) + 6.0 * ((PA * s.B).sum(1)) + 3.0 * ((s.P * s.P * s.R).sum(1))
            if o["k4"] == "path":
                Y += (PA * s.w2[None, :]) @ s.A.T; flops += 2 * n ** 3
            if o["hub"] and not last:
                LA = 2.0 * PA * s.w2[None, :] + 2.0 * s.P * s.B
                LP = (s.A * s.A) * s.w2[None, :] + 2.0 * s.A * s.B + 2.0 * s.P * s.R
                LB = 2.0 * PA
                LR = s.P * s.P
                D21 += LA @ s.A.T + LP @ s.P.T + LB @ s.B.T + LR @ s.R.T; flops += 8 * n ** 3
        d = hermite_relu(alpha, o["K"]); D3 = D3 * o["d3scale"]
        K4 = np.zeros(n)
        if o["k4"] == "path" and srcs:
            Creg = C + o["eps"] * np.mean(var) * np.eye(n)
            Z = np.linalg.solve(Creg, Y.T); flops += (8.0 / 3.0) * n ** 3
            K4 = 12.0 * np.einsum("ij,ji->i", Y, Z); K4 -= np.mean(K4)
        if o["k4diag"] is not None:
            K4 = K4 + o["k4diag"][l]
        m_new = sigma * d[0] - D3 * alpha * ph / (6 * var) + K4 * (alpha ** 2 - 1) * ph / (24 * sigma ** 3)
        out[l] = m_new
        if record is not None:
            record[l] = dict(mu=mu, var=var, alpha=alpha, D3=D3.copy(), D21=D21.copy(), m=m_new.copy(), C=C)
        if last:
            break
        w2 = ph / sigma; w3 = -alpha * ph / var; mg = sigma * d[0]
        second = var * ((1 + alpha ** 2) * Phi + alpha * ph) + D3 * ph / (3 * sigma)
        if o["k4var"]:
            second = second - K4 * alpha * ph / (12 * var)
        rho = C / np.outer(sigma, sigma); np.fill_diagonal(rho, 0.0)
        Kh_new = np.zeros_like(C); rk = np.ones_like(rho); fact = 1.0
        for k in range(1, o["K"] + 1):
            rk = rk * rho; fact *= k; Kh_new += np.outer(d[k], d[k]) * rk / fact
        Kh_new *= np.outer(sigma, sigma)
        if o["hub"]:
            Kh_new += 0.5 * (D21 * w2[:, None] * Phi[None, :] + D21.T * Phi[:, None] * w2[None, :])
        # born fourth-cumulant slices of the previous layer's post-activations would enter here as kappa4 of z; the
        # born slices of THIS layer's h feed the next layer: kept in K22h/K31h (memoryless, one layer) below
        if o["k22born"] and l > 0 and K22z is not None:
            Kh_new += 0.25 * K22z * np.outer(w2, w2)
        if o["k31born"] and l > 0 and K31z is not None:
            Kh_new += (K31z * np.outer(w3, Phi) + K31z.T * np.outer(Phi, w3)) / 6.0
        np.fill_diagonal(Kh_new, second - m_new ** 2)
        if o["pert_off"] or o["pert_var"]:
            dg = np.diag(Kh_new).copy(); Kh_new *= (1 + o["pert_off"]); np.fill_diagonal(Kh_new, dg * (1 + o["pert_var"]))
        # (2,2) and (3,1) slices of this layer's post-activation kappa4, first order in the correlations, transported
        # to the next layer's pre-activation by the uncrossed pattern (W.W) K (W.W)^T and (W.W.W) K31 W^T
        e_vec = 2 * mg * (1 - Phi); vh = second - m_new ** 2
        K22h = np.outer(e_vec, e_vec) * C; np.fill_diagonal(K22h, 0.0)
        K31h = (3 * (1 - Phi) * (vh - mg ** 2))[:, None] * Phi[None, :] * C; np.fill_diagonal(K31h, 0.0)
        Wn = W[l + 1].astype(np.float64); W2 = Wn * Wn
        K22z = W2 @ K22h @ W2.T if o["k22born"] else None
        K31z = (W2 * Wn) @ K31h @ Wn.T if o["k31born"] else None
        # ---- exact gates for the transported sources' slices and diagonals (aggregated residual, post-activation space)
        c_sl = Phi * (1 - Phi) - mg * w2                       # slice gate minus Phi^2
        G3 = Phi - mg * w2 - 0.5 * mg ** 2 * alpha * ph / var  # diagonal gate
        Rg = np.zeros((n, n))
        if o["gate_res"] and srcs:
            Rg = (c_sl[:, None] * D21 * Phi[None, :])          # entries (a, c): c_a Phi_c D21_z[a, c]
            np.fill_diagonal(Rg, (G3 - Phi ** 3) * D3 / 3.0)
        for s in srcs:
            s.P *= Phi[:, None]; s.A *= Phi[:, None]; s.B *= Phi[:, None]; s.R *= Phi[:, None]
        if o["window"]:
            srcs = [s for s in srcs if l + 1 - s.born <= o["window"]]
        # ---- newborn
        A = Phi[:, None] * C
        if o["gc"] and o["hub"]:
            B = (Phi[:, None] * D21.T) * w3[None, :] + 0.5 * (w2[:, None] * D21) * w2[None, :]
        else:
            B = np.zeros((n, n))
        # exact (2,1) slice S21[a, c] = sum_k (1/k!) g_k(a) dz_k(c) C_ac^k,  g_k = F_k(a) - 2 m_a dz_k(a)
        KS = o["KS"]; F = relu2_hermite(alpha, sigma, KS); dz = relu_hermite_z(alpha, sigma, KS)
        g = F - 2 * mg[None, :] * dz
        S21 = np.zeros((n, n)); Ck = np.ones((n, n)); fact = 1.0
        Coff = C.copy(); np.fill_diagonal(Coff, 0.0)
        for k in range(1, KS + 1):
            Ck = Ck * Coff; fact *= k; S21 += np.outer(g[k], dz[k]) * Ck / fact
        # subtract the star's and GC's own (a, a, c) entries: star 2 w2_a A_aa A_ac + w2_c A_ac^2 ; GC (S_3 of e x a x b):
        #   j = a: 2 A_aa B_ac + 2 B_aa A_ac ;  j = c: 2 A_ac B_ac   (entries with the centre at a or at c)
        Aaa = np.diag(A); Baa = np.diag(B)
        own = 2.0 * (w2 * Aaa)[:, None] * A + w2[None, :] * A ** 2
        if o["gc"]:
            own += 2.0 * Aaa[:, None] * B.T + 2.0 * Baa[:, None] * A + 2.0 * A * B.T
        R = (S21 - own) if o["res"] else np.zeros((n, n))
        k3d = sigma ** 3 * k3_exact(alpha)
        np.fill_diagonal(R, (k3d - 3.0 * w2 * Aaa ** 2 - 6.0 * Aaa * Baa) / 3.0)
        R = R.T + Rg.T        # leg columns index the centre a: R_leg[:, a] = residual row a
        srcs.append(Src(np.eye(n), A, B, R, w2, l))
        m, Kh = m_new, Kh_new
    return out, flops
