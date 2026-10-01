"""T17: second trace channel t^mu_l = E[(mu^.a~)^2 a~], mu^ = mu_l/|mu_l|, and the Bolthausen-corrected transport
  t_{l+1} ?= src + 1/2 s2 (sum gamma) Phi*(W^T t_l) + 1/2 c_mu Phi*(W^T t^mu_l) + B phi/s,
  c_mu = sum_a gamma_a m_a^2 / |mu_l|^2 - s2 sum_a gamma_a  (gates depend on w_a through the revealed m_a = mu.w_a).
Exact-centred two-pass MC (N=131072) on bench MLP 0 at n = 1024; regress MC t_{l+1} with and without the new column."""
import sys, numpy as np
from common import bench, phi, Phi
name, i, N = "w1024_d16", 0, int(sys.argv[1]) if len(sys.argv) > 1 else 131072; chunk = 8192
S = bench.load_set(name); W = bench.weights(S, i); L, n, _ = W.shape; s2 = 2.0 / n
def passes(fn):
    rng = np.random.default_rng(99)
    for c in range(N // chunk):
        h = rng.standard_normal((chunk, n)).astype(np.float32)
        for l in range(L):
            z = h @ W[l]; h = np.maximum(z, 0.0); fn(l, z, h)
sa = np.zeros((L, n)); sz = np.zeros((L, n)); sz2 = np.zeros((L, n))
def p1(l, z, h): sa[l] += h.sum(0, dtype=np.float64); sz[l] += z.sum(0, dtype=np.float64); sz2[l] += (z.astype(np.float64) ** 2).sum(0)
passes(p1); mu = sa / N; mz = sz / N; vz = sz2 / N - mz ** 2
muh = mu / np.linalg.norm(mu, axis=1, keepdims=True)
t = np.zeros((L, n)); tm = np.zeros((L, n))
def p2(l, z, h):
    u = h.astype(np.float64) - mu[l]; q = (u ** 2).sum(1); r = (u @ muh[l]) ** 2
    t[l] += (q[:, None] * u).sum(0); tm[l] += (r[:, None] * u).sum(0)
passes(p2); t /= N; tm /= N
np.savez("results/t17_w1024_0.npz", t=t, tm=tm, mu=mu, mz=mz, vz=vz)
Wd = W.astype(np.float64)
import tc
from tc import _xg, _wg, _He
_, states, _ts = tc.predict(W, ret_C="states")
He = _He(6); fact = np.cumprod([1.0] + list(range(1, 7)))
def src_of(m_, S_):
    s_ = np.sqrt(np.diag(S_)); a_ = m_ / s_; R = S_ / np.outer(s_, s_); muG_ = m_ * Phi(a_) + s_ * phi(a_)
    zq = s_[:, None] * (a_[:, None] + _xg[None, :]); r = np.maximum(zq, 0); f = (r - muG_[:, None]) ** 2
    Fk = np.stack([(f * He[k][None]) @ _wg for k in range(7)], 1); Gk = np.stack([(r * He[k][None]) @ _wg for k in range(7)], 1)
    Ro = R.copy(); np.fill_diagonal(Ro, 0.0); Rk = np.ones_like(R); out = np.zeros(len(m_))
    for k in range(1, 7): Rk = Rk * Ro; out += (Rk.T @ Fk[:, k]) * Gk[:, k] / fact[k]
    return out + ((f - Fk[:, :1]) * (r - Gk[:, :1])) @ _wg
for l in range(1, L):
    m = mz[l]; s = np.sqrt(vz[l]); a = m / s; P = Phi(a); p = phi(a); muG = mu[l]
    gam = 2 * P - 2 * muG * p / s; y = Wd[l].T @ t[l - 1]; ym = Wd[l].T @ tm[l - 1]
    cmu = (gam * m ** 2).sum() / (mu[l - 1] @ mu[l - 1]) - s2 * gam.sum()
    yt = t[l]
    # src unknown here: regress on columns, compare R^2 with / without t^mu column (src absorbed via a free column set)
    base = [src_of(*states[l]), P * y, p / s]
    X0 = np.stack(base, 1); X1 = np.stack(base + [P * ym], 1)
    r0 = yt - X0 @ np.linalg.lstsq(X0, yt, rcond=None)[0]; c1 = np.linalg.lstsq(X1, yt, rcond=None)[0]; r1 = yt - X1 @ c1
    A_th = 0.5 * s2 * gam.sum()
    print(f"layer {l+1}: |t^mu|/|t| {np.linalg.norm(tm[l-1])/np.linalg.norm(t[l-1]):.3f}  resid var w/o {r0.var()/yt.var():.4f} with {r1.var()/yt.var():.4f}  "
          f"coef src {c1[0]:.2f} t: {c1[1]:.3f} (th {A_th:.3f})  coef t^mu: {c1[3]:.3f} (th {0.5*cmu:.3f})", flush=True)
