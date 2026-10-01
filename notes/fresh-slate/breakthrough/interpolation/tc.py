"""Trace-channel closure (TC): Gaussian closure (Hermite order 2 covariance) + one n-vector per layer,
t_l = E[|a~_l|^2 a~_l] = Cov(|a~_l|^2, a_l), the chaos-1 (in the weights) part of all third cumulants.
Per layer (z = a_{l-1} W):
  y = W^T t_{l-1};  kappa3(z_p) ~ 3 sigma^2 y_p (diag);  kappa3(z_a, z_a, z_b) ~ sigma^2 y_b (rank-one (2,1) slice)
  first-order Edgeworth injection into E a, E a^2 and Cov(a)  (rank-2 update);
  t_l = src_G(R) + (sigma^2/2) sum(gamma) Phi*y + (sigma^2/2) (t_{l-1} . W u) phi/s,  gamma = 2Phi - 2 mu phi/s, u = 2 mu (1-Phi).
Extra cost over the closure: O(K n^2) per layer."""
import numpy as np
from common import phi, Phi, chi_mean_ratio
from t2_closures import hk

_xg, _wg = np.polynomial.hermite_e.hermegauss(48); _wg = _wg / _wg.sum()
def _He(K):
    H = [np.ones_like(_xg), _xg.copy()]
    for k in range(2, K + 1): H.append(_xg * H[-1] - (k - 1) * H[-2])
    return H

def predict(W, Ksrc=6, inject=True, cov_inject=True, chi=False, Bterm=True, ret_state=False, ret_feats=False, X4=None, ret_C=False, mu_channel=False, x_channel=False, x_scale=1.0):
    W = W.astype(np.float64); L, n, _ = W.shape; sig2 = 2.0 / n
    He = _He(Ksrc); fact = np.cumprod([1.0] + list(range(1, Ksrc + 1)))
    out = []; t = np.zeros(n); tmu = np.zeros(n); mu = None; C = None; ts = []; Cs = []; states = []; Xl = 0.0; Xs = []
    for l in range(L):
        Wl = W[l]
        if l == 0: m = np.zeros(n); S = Wl.T @ Wl; y = np.zeros(n)
        else: m = mu @ Wl; S = Wl.T @ C @ Wl; y = Wl.T @ t; ym = Wl.T @ tmu if mu_channel else None
        v = np.maximum(np.diag(S), 1e-300); s = np.sqrt(v); a = m / s; P = Phi(a); p = phi(a)
        muG = m * P + s * p; secG = (m * m + v) * P + m * s * p
        k3 = 3 * sig2 * y if inject else np.zeros(n)
        mu_new = muG - k3 * a * p / (6 * v)
        sec = secG + k3 * p / (3 * s)
        if x_channel: X4 = [x_scale * Xl] * (l + 1)
        k4 = 3 * sig2 ** 2 * X4[l - 1] if (X4 is not None and l > 0) else 0.0
        if l > 0 and X4 is not None:
            mu_new = mu_new + k4 * p * (a * a - 1) / (24 * v * s)
            sec = sec - k4 * a * p / (12 * v)
        R = S / np.outer(s, s); H = hk(a, 2)
        Cn = np.outer(s * H[0], s * H[0]) * R + 0.5 * np.outer(s * H[1], s * H[1]) * R * R
        if inject and cov_inject:
            u1 = p / s; Cn += 0.5 * sig2 * (np.outer(u1, P * y) + np.outer(P * y, u1))
        if l > 0 and X4 is not None:
            u4 = (p / s) * sig2 * np.sqrt(max(X4[l - 1], 0.0)) / 2; Cn += np.sign(X4[l - 1]) * np.outer(u4, u4)
        np.fill_diagonal(Cn, np.maximum(sec - mu_new * mu_new, 1e-12))
        Cs.append(Cn); states.append((m, S))
        # trace channel of a_l
        zq = s[:, None] * (a[:, None] + _xg[None, :]); r = np.maximum(zq, 0)
        f = (r - muG[:, None]) ** 2
        Fk = np.stack([(f * He[k][None]) @ _wg for k in range(Ksrc + 1)], 1)
        Gk = np.stack([(r * He[k][None]) @ _wg for k in range(Ksrc + 1)], 1)
        Ro = R.copy(); np.fill_diagonal(Ro, 0.0); Rk = np.ones_like(R); src = np.zeros(n)
        for k in range(1, Ksrc + 1):
            Rk = Rk * Ro; src += (Rk.T @ Fk[:, k]) * Gk[:, k] / fact[k]
        src += ((f - Fk[:, :1]) * (r - Gk[:, :1])) @ _wg
        gam = 2 * P - 2 * muG * p / s
        tn = src + 0.5 * sig2 * gam.sum() * P * y
        if Bterm and l > 0:
            u = 2 * muG * (1 - P); tn += 0.5 * sig2 * (t @ (Wl @ u)) * p / s
        if x_channel:
            # Var|a~|^2 by Mehler (Gaussian part) + chaos-0 kappa4 + chaos-1 kappa3 corrections; X = Var - 2||C||^2
            VG = (Ro_sq := None) or 0.0
            Rk2 = np.ones_like(R); Ro2 = R.copy(); np.fill_diagonal(Ro2, 0.0)
            for k in range(1, Ksrc + 1):
                Rk2 = Rk2 * Ro2; VG += Fk[:, k] @ (Rk2 @ Fk[:, k]) / fact[k]
            VG += (((f - Fk[:, :1]) ** 2) @ _wg).sum()
            uu = 2 * muG * (1 - P)
            Vn = VG + (0.25 * sig2 ** 2 * gam.sum() ** 2 * Xl + sig2 * gam.sum() * (uu @ y) if l > 0 else 0.0)
            Xl = Vn - 2 * (Cn ** 2).sum(); Xs.append(Xl)
        if mu_channel:
            muh = mu_new / np.linalg.norm(mu_new); g = p / s; Sv = S @ (P * muh)
            tmn = 2 * P * (S @ (muh * g * Sv)) + g * Sv ** 2
            if l > 0:
                muh_prev = mu / np.linalg.norm(mu)
                cmu = (gam * m * m).sum() / (mu @ mu) - sig2 * gam.sum()
                tn = tn + 0.5 * cmu * P * ym
                nu = Wl @ (P * muh); c = nu @ muh_prev; nperp2 = nu @ nu - c * c
                tmn = tmn + P * (c * c * ym + (nperp2 / n) * y)
            tmu = tmn
        t = tn; mu = mu_new; C = Cn; out.append(mu); ts.append(t)
    out = np.stack(out)
    if chi: out = out * chi_mean_ratio(n)
    if ret_C == 'X': return out, np.array(Xs)
    if ret_C == 'states': return out, states, np.stack(ts)
    if ret_C: return out, Cs
    if ret_feats: return out, a, s
    return (out, np.stack(ts)) if ret_state else out
