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


def predict(W, mode="lin", R=6, use_k4=True, mean_var="edge", d21=False, win=1):
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
                for s in range(l if mode == "lin" else max(l - win, 0)):
                    K = Px[s].T @ b                          # Cov of the linear legs, (k, j)
                    M = Pfrom[s] * cs[s][:, None] * K
                    k3 += 6 * (M * K).sum(0)
                    if d21:   # kappa(y_j, y_j, y_m) = 2 b_j^T Q_m b_j + 4 b_j^T Q_j b_m
                        D21 += 2 * (K * K).T @ (Pfrom[s] * cs[s][:, None]) + 4 * M.T @ K
                    if use_k4:
                        V += Px[s] @ M
                k4 = 48 * (V * V).sum(0) if use_k4 else np.zeros(n)
                if mode == "hyb":    # the win most recent births on renormalised legs (closure covariance of the birth layer)
                    for s in range(max(l - win, 0), l):
                        Pc = Pfrom[s] * cs[s][:, None]
                        K = (Cs[s] * betas[s][None, :]) @ Pfrom[s]
                        M = Pc * K
                        k3 = k3 + 6 * (M * K).sum(0)
                        if d21:
                            D21 = D21 + 2 * (K * K).T @ Pc + 4 * M.T @ K
                        if use_k4:
                            k4 = k4 + 48 * (M * (Cs[s] @ M)).sum(0)
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
        Ea2 = second_moment(mu, var, k3, k4) if mean_var == "edge" else (mu * mu + var) * ndtr(t) + mu * s_ * _phi(t)
        np.fill_diagonal(Ca, np.maximum(Ea2 - Ea * Ea, 1e-300))
        Cprev = C
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
