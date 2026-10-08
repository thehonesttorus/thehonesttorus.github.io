# Exact check of the per-unit gate symbol at a Gaussian collective state: for e_a = relu(mu + u.t + s*xi) - mean,
# t ~ N(0, Lam), xi ~ N(0,1) independent, the symbol N_a = kappa(e_a, e_a, t, t) is rank one along Lam u:
#     N_a = c_a (Lam u)(Lam u)^T,   c_a = E[(f^2)''] - 2 E[f']^2 - 2 E[f] E[f''] + E[g''],
# with f(p) = E[relu(mu + p + s xi)], g(p) = Var(relu(mu + p + s xi)), p = u.t ~ N(0, u^T Lam u). Gauss-Hermite in (t, xi).
import numpy as np
from math import erf, sqrt, pi, exp
rng = np.random.default_rng(5)
K = 2; mu, s = 0.3, 0.8
A = rng.normal(size=(K, K)); Lam = A @ A.T + 0.2 * np.eye(K); u = rng.normal(size=K)
x, w = np.polynomial.hermite_e.hermegauss(80); w = w / w.sum()
Lc = np.linalg.cholesky(Lam)
# tensor grid over (g1, g2, xi): t = Lc g
G1, G2, XI = np.meshgrid(x, x, x, indexing="ij"); WW = np.einsum("i,j,k->ijk", w, w, w)
T = np.einsum("km,mijl->kijl", Lc, np.stack([G1, G2]))
y = np.maximum(mu + np.einsum("k,kijl->ijl", u, T) + s * XI, 0.0)
E = lambda F: float(np.sum(WW * F))
e = y - E(y)
def cum(a, b, c, d):
    return E(a * b * c * d) - E(a * b) * E(c * d) - E(a * c) * E(b * d) - E(a * d) * E(b * c)
N = np.array([[cum(e, e, T[k], T[m]) for m in range(K)] for k in range(K)])
v = Lam @ u
# f, g as functions of p (1-D quadrature for the derivatives' expectations)
Phi = lambda z: 0.5 * (1 + erf(z / sqrt(2))); phi = lambda z: exp(-z * z / 2) / sqrt(2 * pi)
def fm(p):
    m = mu + p; a = m / s; return m * Phi(a) + s * phi(a)
def gv(p):
    m = mu + p; a = m / s; m1 = m * Phi(a) + s * phi(a); m2 = (m * m + s * s) * Phi(a) + m * s * phi(a); return m2 - m1 * m1
sp = sqrt(u @ Lam @ u); h = 1e-3
pp = sp * x
d1 = lambda F, p: (F(p + h) - F(p - h)) / (2 * h); d2 = lambda F, p: (F(p + h) - 2 * F(p) + F(p - h)) / h ** 2
Ef1 = sum(wi * d1(fm, p) for wi, p in zip(w, pp)); Ef = sum(wi * fm(p) for wi, p in zip(w, pp))
Ef2 = sum(wi * d2(fm, p) for wi, p in zip(w, pp)); Eff2 = sum(wi * d2(lambda q: fm(q) ** 2, p) for wi, p in zip(w, pp))
Eg2 = sum(wi * d2(gv, p) for wi, p in zip(w, pp))
c = Eff2 - 2 * Ef1 ** 2 - 2 * Ef * Ef2 + Eg2
print("N_a (quadrature)\n", N, "\nc_a (Lam u)(Lam u)^T\n", c * np.outer(v, v))
print("rank-one residual / |N|: %.2e" % (np.linalg.norm(N - c * np.outer(v, v)) / np.linalg.norm(N)))
