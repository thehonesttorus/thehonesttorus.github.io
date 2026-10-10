"""The collective-conditional closure (CC1): Stage 18/21/23's collective form as an analytic forward state.

At every layer the post-activation law is represented as a mixture over the collective coordinate c = mu_hat . h of
conditional Gaussians (Stage 22 audit, Prop. 5.1: the attribution uses the actual conditional mean and variance):
    h | c ~ N( mu + b (c - cbar),  rho(c) K_t ),   c ~ N(cbar, V_c),
    b = K mu_hat / V_c                      (regression of h on c; b . mu_hat = 1)
    K_t = K - V_c b b^T                     (law of total covariance; K_t mu_hat = 0)
    rho(c) = 1 + gamma (c - cbar) / sd_c    (heteroscedastic transverse scale; gamma = Cov(c, r^2) / (sd_c E r^2)),
with r^2 = |P_perp h|^2. The plain Gaussian closure is gamma = 0 with the mixture collapsed; it cannot carry the
collective scale mode, which relu passes multiplicatively (homogeneity), and it under-reports Var c by 26-45%
(notes/stage21/results/tier2_oracle.txt). Here each Gauss-Hermite node of c is pushed through the layer exactly
within the model: pre-activation z | c ~ N(W mu + (W b) dc, rho(c) W K_t W^T), relu moments by the Price/Hermite
series with the shared correlation of W K_t W^T and node-specific standardised means, and the next layer's
(mu, K, V_c, Cov(c, r^2)) by the law of total (co)variance over nodes. Collective second and third moments come
from Hermite cross series of relu and relu^2 (exact for a bivariate Gaussian pair), so nothing is fitted.

Cost per layer (float32 under flopscope, notes/compute/COST_MODEL.md): one covariance sandwich (3.33 n^3 with a
Cholesky factor) plus about 3 Q K n^2 for Q nodes and K Hermite terms (Q = 7, K = 10: 0.2 n^3).
Convention: W[l] of shape (n_out, n_in), z = W h.
"""
import numpy as np
from scipy.special import ndtr

SQ2PI = np.sqrt(2 * np.pi)


def _phi(a): return np.exp(-0.5 * a * a) / SQ2PI


def trunc_moments(a, pmax=4):
    """M_p(a) = E (u + a)_+^p, u ~ N(0, 1), p = 0..pmax (M_p = a M_(p-1) + (p-1) M_(p-2))."""
    P, ph = ndtr(a), _phi(a); M = [P, ph + a * P]
    for p in range(2, pmax + 1): M.append(a * M[p - 1] + (p - 1) * M[p - 2])
    return M


def hermite_coeffs(a, K):
    """Hermite coefficients E[f He_t] of f = (u + a)_+ (d_t) and f = (u + a)_+^2 (e_t), t = 0..K."""
    M = trunc_moments(a, 2); ph, P = _phi(a), ndtr(a); x = -a
    He = [np.ones_like(a), x.copy()]
    for j in range(2, K): He.append(x * He[-1] - (j - 1) * He[-2])
    d = [M[1], P] + [ph * He[t - 2] for t in range(2, K + 1)]
    e = [M[2], 2 * M[1], 2 * P] + [2 * ph * He[t - 3] for t in range(3, K + 1)]
    return d, e


class CC1:
    def __init__(self, Q=7, K=10, collective=True):
        self.Q, self.K, self.collective = Q, K, collective
        g, w = np.polynomial.hermite_e.hermegauss(Q); self.g, self.w = g, w / w.sum()
        self.fact = [float(np.prod(np.arange(1, t + 1))) for t in range(K + 1)]

    def layer(self, W, st):
        """One layer. st: dict(mu, K, Vc, gamma) of the incoming post-activation law (st=None: input N(0, I))."""
        K_ = self.K
        if st is None:                                       # exact Gaussian input: a single node
            mu_z = np.zeros(W.shape[0]); Ct = W @ W.T; nodes = [(1.0, 0.0, 1.0)]; Wb = np.zeros(W.shape[0])
        else:
            mu, Kh, Vc = st["mu"], st["K"], st["Vc"]
            muh = mu / np.linalg.norm(mu); Km = Kh @ muh; b = Km / Vc
            mu_z = W @ mu; Wb = W @ b
            Cz = W @ Kh @ W.T; Ct = Cz - Vc * np.outer(Wb, Wb)       # W K_t W^T as a rank-one downdate
            sdc = np.sqrt(Vc); trKt = np.trace(Kh) - Vc * (b @ b)       # tr K_t = E tr Cov(h | c)
            gam = st["cov_cr2"] / (sdc * trKt) if self.collective else 0.0
            if self.collective:
                nodes = [(wq, sdc * gq, max(1.0 + gam * gq, 1e-3)) for gq, wq in zip(self.g, self.w)]
            else:
                nodes = [(1.0, 0.0, 1.0)]; Ct = Cz                     # plain Gaussian closure
        st_sd = np.sqrt(np.clip(np.diag(Ct), 1e-300, None)); R = Ct / np.outer(st_sd, st_sd); np.fill_diagonal(R, 0.0)
        # per node: standardised means and Hermite coefficients
        per = []
        for wq, dc, rq in nodes:
            sd = st_sd * np.sqrt(rq); a = (mu_z + Wb * dc) / sd
            d, e = hermite_coeffs(a, K_); M = trunc_moments(a, 4)
            per.append(dict(w=wq, sd=sd, a=a, d=d, e=e, M=M, m=sd * M[1], s2=sd * sd * M[2]))
        mu_n = sum(p["w"] * p["m"] for p in per); muh_n = mu_n / np.linalg.norm(mu_n)
        # covariance: within-node Hermite series (shared R, rank-Q coefficient sums) + between-node spread
        Kn = np.zeros_like(R); Rk = np.ones_like(R)
        cross = {"Vc": [0.0] * len(per), "Ccq": [0.0] * len(per)}
        for t in range(1, K_ + 1):
            Rk = Rk * R
            A = sum(p["w"] * np.outer(p["sd"] * p["d"][t], p["sd"] * p["d"][t]) for p in per)
            Kn += Rk * A / self.fact[t]
            for i, p in enumerate(per):                      # node quadratic forms for the collective moments
                u = muh_n * p["sd"] * p["d"][t]; v = p["sd"] ** 2 * p["e"][t]
                Ru = Rk @ u
                cross["Vc"][i] += (u @ Ru) / self.fact[t]
                cross["Ccq"][i] += (v @ Ru) / self.fact[t]
        # exact diagonals: Var(relu z) and Cov(relu z, relu^2 z) per node
        diag_w = np.zeros_like(mu_n)
        for i, p in enumerate(per):
            M = p["M"]; sd = p["sd"]
            var_d = sd ** 2 * (M[2] - M[1] ** 2); cov_d = sd ** 3 * (M[3] - M[1] * M[2])
            diag_w += p["w"] * var_d
            cross["Vc"][i] += np.sum(muh_n ** 2 * var_d); cross["Ccq"][i] += np.sum(muh_n * cov_d)
        np.fill_diagonal(Kn, diag_w)
        for p in per:                                        # between-node covariance
            dm = p["m"] - mu_n; Kn += p["w"] * np.outer(dm, dm)
        # collective moments of the new layer: c' = muh_n . h', q' = |h'|^2, r'^2 = q' - c'^2
        cbar = muh_n @ mu_n
        Ec = np.array([muh_n @ p["m"] for p in per]); Vcq = np.array(cross["Vc"]); Ccq = np.array(cross["Ccq"])
        Eq = np.array([np.sum(p["s2"]) for p in per]); ws = np.array([p["w"] for p in per])
        Er2_q = Eq - Ec ** 2 - Vcq                            # E[r'^2 | node]
        Ccr2_q = Ccq - 2 * Ec * Vcq                           # Cov(c', r'^2 | node), delta method on c'^2
        Vc_n = muh_n @ Kn @ muh_n
        Er2 = ws @ Er2_q
        cov_cr2 = ws @ (Ccr2_q + (Ec - cbar) * (Er2_q - Er2))
        pre = dict(mu_z=mu_z, Wb=Wb, vt=np.diag(Ct).copy(), nodes=nodes)    # this layer's pre-activation mixture
        return dict(mu=mu_n, K=Kn, Vc=Vc_n, cbar=cbar, Er2=Er2, cov_cr2=cov_cr2, pre=pre)

    def run(self, Ws):
        st, out, coll = None, [], []
        for Wl in Ws:
            st = self.layer(Wl, st); out.append(st["mu"])
            coll.append((st["cbar"], st["Vc"], st["cov_cr2"], st["Er2"]))
        self.last = st
        return np.array(out), coll
