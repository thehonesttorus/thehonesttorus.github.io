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

def predict(W, Ksrc=6, inject=True, cov_inject=True, chi=False, Bterm=True, ret_state=False):
    W = W.astype(np.float64); L, n, _ = W.shape; sig2 = 2.0 / n
    He = _He(Ksrc); fact = np.cumprod([1.0] + list(range(1, Ksrc + 1)))
    out = []; t = np.zeros(n); mu = None; C = None; ts = []
    for l in range(L):
        Wl = W[l]
        if l == 0: m = np.zeros(n); S = Wl.T @ Wl; y = np.zeros(n)
        else: m = mu @ Wl; S = Wl.T @ C @ Wl; y = Wl.T @ t
        v = np.maximum(np.diag(S), 1e-300); s = np.sqrt(v); a = m / s; P = Phi(a); p = phi(a)
        muG = m * P + s * p; secG = (m * m + v) * P + m * s * p
        k3 = 3 * sig2 * y if inject else np.zeros(n)
        mu_new = muG - k3 * a * p / (6 * v)
        sec = secG + k3 * p / (3 * s)
        R = S / np.outer(s, s); H = hk(a, 2)
        Cn = np.outer(s * H[0], s * H[0]) * R + 0.5 * np.outer(s * H[1], s * H[1]) * R * R
        if inject and cov_inject:
            u1 = p / s; Cn += 0.5 * sig2 * (np.outer(u1, P * y) + np.outer(P * y, u1))
        np.fill_diagonal(Cn, np.maximum(sec - mu_new * mu_new, 1e-12))
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
        t = tn; mu = mu_new; C = Cn; out.append(mu); ts.append(t)
    out = np.stack(out)
    if chi: out = out * chi_mean_ratio(n)
    return (out, np.stack(ts)) if ret_state else out
