"""Identity checks for notes/stage10/s10_kikuchi.tex (Kikuchi sparsification, dualities, two-tier dynamics). No estimator runs.
  python scripts/verify_kikuchi.py"""
import numpy as np
from math import sqrt, pi, acos, asin
from scipy.integrate import quad
from scipy.stats import norm, multivariate_normal
rng = np.random.default_rng(12)
J = lambda th: np.sin(th) + (pi - th) * np.cos(th)

print("1. Schur-Weyl reading: complete-graph Kikuchi Laplacian on E_k has eigenvalue k(n+1-k) = C(n,2) - (content sum of (n-k,k))")
for n in (8, 12):
    for k in range(0, n // 2 + 1):
        content = sum(range(n - k)) + sum(j - 1 for j in range(k))
        assert n * (n - 1) // 2 - content == k * (n + 1 - k)
print("   holds for n = 8, 12 and all k")

print("2. gate covariance = integrated wall-current Gram along Eldan's path (layer 1, any means): Cov(g_a,g_b) vs integral")
for (ma, mb, rho) in ((0.0, 0.0, 0.3), (0.7, -0.4, 0.6), (1.5, 1.2, -0.5), (0.8, 0.8, 1.0)):
    # z_a = m_a + s_a, unit variances, corr rho. Along Eldan's path the localized law of (z_a, z_b) at correlation-time s = t/(1+t):
    # integrated current Gram = int_0^1 d/ds P(z_a>0, z_b'>0) where (z_a, z_b') correlated rho*s  (Price's identity)
    def P2(r):
        if abs(r) >= 1: return norm.cdf(min(ma, mb)) if r > 0 else max(0.0, norm.cdf(ma) + norm.cdf(mb) - 1)
        return multivariate_normal([0, 0], [[1, r], [r, 1]]).cdf([ma, mb])
    cov = P2(rho) - norm.cdf(ma) * norm.cdf(mb)
    gram = quad(lambda s: rho * np.exp(-(ma ** 2 + mb ** 2 - 2 * rho * s * ma * mb) / (2 * (1 - rho ** 2 * s ** 2) + 1e-300)) / (2 * pi * np.sqrt(max(1 - rho ** 2 * s ** 2, 1e-300))), 0, 1, limit=200)[0]
    print("   m=(%.1f,%.1f) rho=%.1f: Cov %.6f   integrated current Gram %.6f" % (ma, mb, rho, cov, gram))

print("3. trickle in depth: fraction of layer-l gates on which two conditionally independent copies disagree = theta_{l-1}/pi")
n, L, P = 512, 16, 64
Ws = [rng.standard_normal((n, n)) * sqrt(2 / n) for _ in range(L)]
for rho in (0.5, 0.9, 0.99):
    th = [acos(rho)]
    for l in range(1, L): th.append(acos(min(1.0, J(th[-1]) / pi)))
    X = rng.standard_normal((P, n)); Xp = rho * X + sqrt(1 - rho ** 2) * rng.standard_normal((P, n))
    h, hp, dis = X, Xp, []
    for W in Ws:
        z, zp = h @ W.T, hp @ W.T; dis.append(np.mean((z > 0) != (zp > 0))); h, hp = np.maximum(z, 0), np.maximum(zp, 0)
    rows = " ".join("%d:%.3f/%.3f" % (l + 1, dis[l], th[l] / pi) for l in (0, 1, 3, 7, 11, 15))
    print("   rho=%.2f (t=%.1f)  layer: measured/predicted  %s" % (rho, rho / (1 - rho), rows))
