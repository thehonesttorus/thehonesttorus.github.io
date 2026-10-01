"""Facet-birth transport (FBT), faces design v1: numpy prototype for Stage Q.

State per layer l (pre-activation z_l):  mu (n), C (n x n), plus the face-averaged arrows.
  1. Face measures of each neuron from (mu, sigma): mass Phi(t), facet density phi(t)/sigma and its derivatives
     (E2: the Hermite profile of ReLU is the sequence of face measures of increasing codimension).
  2. Covariance of the activations by the face-measure Mehler series
     Cov(a_i, a_k) = sum_{r=1..R} fh_i(r) fh_k(r) rho_ik^r / r!,  fh(1) = sigma Phi(t), fh(r) = sigma He_{r-2}(-t) phi(t);
     exact diagonal.  Arrow: mu' = W^T E a, C' = W^T Cov(a) W.
  3. Non-Gaussianity: facet births on the Gaussian legs of the input, transported by the face-averaged arrows (E5):
     z_l ~= x b_l + x^T Q_l x,  b_l = P_{x->l}[:, j],  Q_{l,j} = sum_{s<l} sum_k P_{s->l}[k,j] c_{s,k} v_{s,k} v_{s,k}^T,
     v_{s,k} = P_{x->s}[:, k],  c = facet density / 2 = phi(t) / (2 sigma),  beta = Phi(t).
     kappa3 = 6 b^T Q b,  kappa4 = 48 b^T Q^2 b  (loop terms tr Q^3, tr Q^4 dropped: n^{-1/2}-suppressed).
     'mem' variant: only s = l-1, legs from the closure covariance (the one-step Wick tree).
  4. Readout: Edgeworth (facet-density form, DESIGN.md section 3) with (mu, C_jj, kappa3, kappa4).
"""
import numpy as np
from scipy.special import ndtr, eval_hermitenorm, factorial

SQ2PI = np.sqrt(2 * np.pi)


def _phi(t):
    return np.exp(-0.5 * t * t) / SQ2PI


def readout(mu, var, k3, k4):
    s = np.sqrt(var); t = mu / s; f = _phi(t)
    base = mu * ndtr(t) + s * f
    return base - k3 * t * f / (6 * s ** 2) + k4 * (t * t - 1) * f / (24 * s ** 3) \
        + k3 ** 2 * (t ** 4 - 6 * t * t + 3) * f / (72 * s ** 5)


def second_moment(mu, var, k3, k4):
    """E relu(z)^2 with Edgeworth corrections: int (s u + mu)_+^2 He_k phi = 2 s^2 int 1{u>-t} He_{k-2} phi = 2 s^2 He_{k-3}(-t) phi(t)."""
    s = np.sqrt(var); t = mu / s; f = _phi(t)
    base = (mu * mu + var) * ndtr(t) + mu * s * f
    h0 = 1.0; h1 = -t; h3 = -t ** 3 + 3 * t
    return base + k3 / (6 * s ** 3) * 2 * s * s * h0 * f + k4 / (24 * s ** 4) * 2 * s * s * h1 * f \
        + k3 ** 2 / (72 * s ** 6) * 2 * s * s * h3 * f


def _rmom(n):
    from scipy.special import gammaln
    return [np.exp(gammaln((n + k) / 2) - gammaln(n / 2)) * (2.0 / n) ** (k / 2) for k in range(5)]


def radial_dk4(mu, var, k3, n):
    """Exact radial-mode excess of kappa4 for z = R y, R = |x|/sqrt(n) independent of the sphere-level y
    (positive homogeneity): y taken with the closure's first three cumulants and kappa4(y) = 0."""
    R = _rmom(n)
    def k4_of(r1, r2, r3, r4):
        m1y = mu / r1; m2y = mu * mu + var            # E R^2 = 1 exactly
        vy = m2y - m1y ** 2
        m3y = k3 + 3 * m1y * vy + m1y ** 3
        m4y = 3 * vy ** 2 + 4 * k3 * m1y + 6 * m1y ** 2 * vy + m1y ** 4
        m1, m2, m3, m4 = r1 * m1y, r2 * m2y, r3 * m3y, r4 * m4y
        return m4 - 4 * m3 * m1 - 3 * m2 ** 2 + 12 * m2 * m1 ** 2 - 6 * m1 ** 4
    return k4_of(R[1], R[2], R[3], R[4]) - k4_of(1.0, 1.0, 1.0, 1.0)


def mehler_cov(mu, C, R=6):
    s = np.sqrt(np.maximum(np.diag(C), 1e-300)); t = mu / s; f = _phi(t)
    rho = C / np.outer(s, s)
    out = np.zeros_like(C)
    rp = np.ones_like(C)
    for r in range(1, R + 1):
        fh = s * ndtr(t) if r == 1 else s * eval_hermitenorm(r - 2, -t) * f
        rp = rp * rho
        out += np.outer(fh, fh) * rp / factorial(r)
    return out


def predict(W, mode="lin", R=6, use_k4=True, mean_var="edge", d21=False, win=1, k22=False, winonly=False, radial=False, cpass=0.0):
    """mode: 'gauss' (no kappa3/4), 'lin' (FBT: facet births on Gaussian input legs, all depths), 'mem' (one-step tree)."""
    W = np.asarray(W, dtype=np.float64)
    L, n, _ = W.shape
    out = []
    mu = np.zeros(n); C = W[0].T @ W[0]
    Px = [W[0]]                 # x-space face-averaged Jacobians of the linear part, P_{x->l}
    Pfrom = []                  # Pfrom[s] = P_{s->l} for the current l, s < l
    cs = []                     # facet coefficients per earlier layer
    Cs = []; betas = []         # closure covariances and slopes per earlier layer
    k3 = np.zeros(n); k4 = np.zeros(n)
    for l in range(L):
        if l > 0:
            mu = Ea @ W[l]
            C = W[l].T @ Ca @ W[l]
            Px.append((Px[-1] * beta[None, :]) @ W[l])
            Pfrom = [(P * beta[None, :]) @ W[l] for P in Pfrom] + [W[l]]
            cs.append(cfac); Cs.append(Cprev); betas.append(beta)
            if mode in ("lin", "hyb"):
                b = Px[l]
                k3 = np.zeros(n); V = np.zeros((n, n)); D21 = np.zeros((n, n))
                for s in range(l if mode == "lin" else (0 if winonly else max(l - win, 0))):
                    K = Px[s].T @ b                          # Cov of the linear legs, (k, j)
                    M = Pfrom[s] * cs[s][:, None] * K
                    k3 += 6 * (M * K).sum(0)
                    if d21:   # kappa(y_j, y_j, y_m) = 2 b_j^T Q_m b_j + 4 b_j^T Q_j b_m
                        D21 += 2 * (K * K).T @ (Pfrom[s] * cs[s][:, None]) + 4 * M.T @ K
                    if use_k4:
                        V += Px[s] @ M
                k4 = 48 * (V * V).sum(0) if use_k4 else np.zeros(n)
                if mode == "hyb":    # the win most recent births on renormalised legs (closure covariance of the birth layer)
                    K22 = np.zeros((n, n)); K31 = np.zeros((n, n))
                    for s in range(max(l - win, 0), l):
                        Pc = Pfrom[s] * cs[s][:, None]
                        K = (Cs[s] * betas[s][None, :]) @ Pfrom[s]
                        M = Pc * K
                        k3 = k3 + 6 * (M * K).sum(0)
                        if d21:
                            D21 = D21 + 2 * (K * K).T @ Pc + 4 * M.T @ K
                        if use_k4:
                            k4 = k4 + 48 * (M * (Cs[s] @ M)).sum(0)
                        if k22:   # two-birth chains with a diagonal inner covariance
                            dC = np.diag(Cs[s])[:, None]
                            A = (K * K).T @ (dC * Pc * Pc)                # b_j Q_m^2 b_j
                            B = M.T @ (dC * M)                            # b_j Q_j Q_m b_m
                            K22 += 8 * (A + A.T) + 32 * B
                            K31 += 24 * (M * K).T @ (dC * Pc) + 24 * (M * Pc).T @ (dC * K)
            elif mode == "mem":
                K = (Cprev * beta[None, :]) @ W[l]          # Cov(z_{l-1,k}, linear image at l)
                M = W[l] * cfac[:, None] * K
                k3 = 6 * (M * K).sum(0)
                if d21:
                    D21 = 2 * (K * K).T @ (W[l] * cfac[:, None]) + 4 * M.T @ K
                if use_k4:
                    # b^T Q^2 b with legs in z_{l-1} space: Q_j = sum_k W_kj c_k e_k e_k^T, b = beta W (Gaussian z_{l-1}, cov Cprev)
                    Y = Cprev @ M
                    k4 = 48 * (M * Y).sum(0)
                else:
                    k4 = np.zeros(n)
        var = np.maximum(np.diag(C), 1e-300); s_ = np.sqrt(var); t = mu / s_
        if mode == "gauss":
            k3 = np.zeros(n); k4 = np.zeros(n)
        if cpass and l > 1 and mode == "hyb" and d21:
            # curvature passage of old skew through the facets of layer l-1: 3 sym[w2_i Phi_j w2_k D21(z)_ik C_jk]
            w2 = 2 * cfac
            U = W[l] * w2[:, None]
            V = (w2[:, None] * W[l]) * ((Cprev * beta[None, :]) @ W[l])
            k3 = k3 + cpass * (U * (D21prev @ V)).sum(0)
        if d21 and l > 0 and mode == "hyb":
            D21prev_next = D21
        if radial and l > 0:
            k4 = k4 + radial_dk4(mu, var, k3, n)
        Ea = readout(mu, var, k3, k4)
        out.append(Ea)
        if l == L - 1:
            break
        beta = ndtr(t)
        cfac = _phi(t) / (2 * s_)
        Ca = mehler_cov(mu, C, R)
        if d21 and l > 0 and mode in ("lin", "mem", "hyb"):
            # bivariate Edgeworth, leading order: dE[a_i a_k] = D21_ik phi_i Phi_k / (2 s_i) + D21_ki Phi_i phi_k / (2 s_k)
            A1 = (_phi(t) / s_)[:, None] * beta[None, :] * D21
            Ca = Ca + 0.5 * (A1 + A1.T)
            if k22 and mode == "hyb":
                ph = _phi(t)
                Ca = Ca + 0.25 * K22 * np.outer(ph / s_, ph / s_)
                A3 = (-t * ph / s_ ** 2)[:, None] * beta[None, :] * K31 / 6.0
                Ca = Ca + (A3 + A3.T)
        Ea2 = second_moment(mu, var, k3, k4) if mean_var == "edge" else (mu * mu + var) * ndtr(t) + mu * s_ * _phi(t)
        np.fill_diagonal(Ca, np.maximum(Ea2 - Ea * Ea, 1e-300))
        Cprev = C
        if d21 and l > 0 and mode == "hyb":
            D21prev = D21prev_next
    return np.stack(out)


def predict_gauss(W):
    return predict(W, mode="gauss")


def predict_lin(W):
    return predict(W, mode="lin")


def predict_lin_k3(W):
    return predict(W, mode="lin", use_k4=False)


def predict_mem(W):
    return predict(W, mode="mem")


def predict_mem21(W):
    return predict(W, mode="mem", d21=True)


def predict_lin21(W):
    return predict(W, mode="lin", d21=True)


def predict_hyb21(W):
    return predict(W, mode="hyb", d21=True)


def predict_hyb21_w2(W):
    return predict(W, mode="hyb", d21=True, win=2)


def predict_hyb21_w3(W):
    return predict(W, mode="hyb", d21=True, win=3)


def _hyb(w):
    def f(W):
        return predict(W, mode="hyb", d21=True, win=w)
    f.__name__ = f"predict_hyb21_w{w}"
    return f


for _w in (4, 6, 8, 12, 16):
    globals()[f"predict_hyb21_w{_w}"] = _hyb(_w)


def predict_legs_nod21(W):
    return predict(W, mode="hyb", d21=False, win=16)


def predict_legs_nok4(W):
    return predict(W, mode="hyb", d21=True, win=16, use_k4=False)


def predict_legs_k22(W):
    return predict(W, mode="hyb", d21=True, win=16, k22=True)


def predict_w6_k22(W):
    return predict(W, mode="hyb", d21=True, win=6, k22=True)


def predict_win4(W):
    return predict(W, mode="hyb", d21=True, win=4, winonly=True)


def predict_win6(W):
    return predict(W, mode="hyb", d21=True, win=6, winonly=True)


def predict_w16_rad(W):
    return predict(W, mode="hyb", d21=True, win=16, radial=True)


def predict_win6_rad(W):
    return predict(W, mode="hyb", d21=True, win=6, winonly=True, radial=True)


def predict_gauss_rad(W):
    return predict(W, mode="gauss", radial=True)


def predict_w16_cp3(W):
    return predict(W, mode="hyb", d21=True, win=16, cpass=3.0)


def predict_w16_cp1(W):
    return predict(W, mode="hyb", d21=True, win=16, cpass=1.0)
