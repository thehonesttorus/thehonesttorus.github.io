"""T6: one-step recursion for the trace channel t_l = E[|a~_l|^2 a~_l] (MC values from t5 as oracle inputs).
Derived (Edgeworth to first order in kappa3(z), chaos-0 contraction of the transport Gram):
  t_{l+1,c} ~ src_c + A * Phi_c (W^T t_l)_c + B * phi_c/s_c
  src_c   = sum_a Cov_G((relu z_a - mu_a)^2, relu z_c)  (Gaussian reference, Hermite series in rho_ac, O(n^2))
  A_theory = (sigma^2/2) sum_a gamma_a,  gamma_a = E f_a'' = 2 Phi_a - 2 mu_a phi_a / s_a
Regress MC t_{l+1} on the three columns; report R^2 and fitted vs theory coefficients."""
import sys
import numpy as np
from common import bench, phi, Phi
from t2_closures import hk
name, i = sys.argv[1], int(sys.argv[2])
S_ = bench.load_set(name); W = bench.weights(S_, i).astype(np.float64); T = S_["means"][i]; L, n, _ = W.shape
d = np.load(sorted(__import__("glob").glob(f"results/t5_{name}_{i}_N*.npz"))[-1]); tmc = d["t"]
sig2 = 2.0 / n; K = 6
xg, wg = np.polynomial.hermite_e.hermegauss(60); wg = wg / wg.sum()
He = [np.ones_like(xg), xg]
for k in range(2, K + 1): He.append(xg * He[-1] - (k - 1) * He[-2])
fact = np.cumprod([1] + list(range(1, K + 1)))
# closure states (K=2), keep (m, S) per layer
m = np.zeros(n); Sz = W[0].T @ W[0]; states = []
for l in range(L):
    if l > 0: m = mu @ W[l]; Sz = W[l].T @ C @ W[l]
    v = np.diag(Sz); s = np.sqrt(v); a = m / s; P = Phi(a); p = phi(a)
    mu = m * P + s * p; sec = (m * m + v) * P + m * s * p
    R = Sz / np.outer(s, s); H = hk(a, 2)
    C = np.outer(s * H[0], s * H[0]) * R + 0.5 * np.outer(s * H[1], s * H[1]) * R * R
    np.fill_diagonal(C, sec - mu * mu); states.append((m, Sz, s, a, mu, R))
for l in range(1, L):
    m, Sz, s, a, mu, R = states[l]
    zq = s[:, None] * (a[:, None] + xg[None, :])                  # (n, nodes)
    r = np.maximum(zq, 0)
    f = (r - mu[:, None]) ** 2; g = r
    Fk = np.stack([(f * He[k][None]) @ wg for k in range(K + 1)], 1)   # E[f He_k]
    Gk = np.stack([(g * He[k][None]) @ wg for k in range(K + 1)], 1)
    src = np.zeros(n); Rk = np.ones_like(R); Ro = R.copy(); np.fill_diagonal(Ro, 0)
    for k in range(1, K + 1):
        Rk = Rk * Ro; src += (Rk.T @ Fk[:, k]) * Gk[:, k] / fact[k]
    src += ((f - Fk[:, :1]) * (g - Gk[:, :1])) @ wg                 # a = c term (1-D Gaussian third moment)
    Pc = Phi(a); pc = phi(a)
    X1 = Pc * (W[l].T @ tmc[l - 1]); X3 = pc / s
    gam = 2 * Pc - 2 * mu * pc / s; A_th = 0.5 * sig2 * gam.sum()
    y = tmc[l]
    X = np.stack([src, X1, X3], 1); c, *_ = np.linalg.lstsq(X, y, rcond=None)
    r2 = 1 - ((y - X @ c) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    th = src + A_th * X1
    c3 = np.linalg.lstsq(X3[:, None], y - th, rcond=None)[0][0]
    r2th = 1 - ((y - th - c3 * X3) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    r2src = 1 - ((y - src) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    print(f"layer {l+1}: R2 fit {r2:.3f} (src {c[0]:.2f}, A {c[1]:.3f} vs theory {A_th:.3f}, B {c[2]:+.3f}); "
          f"R2 theory(+B) {r2th:.3f}; R2 src only {r2src:.3f}; rms t {np.sqrt((y**2).mean()):.3f}", flush=True)
