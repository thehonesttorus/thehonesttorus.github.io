# Level-2 sketch: Bayes error E Var(F | w.m, q=u^T Sigma' u)  vs  Gaussian-closure error |E relu(N(w.m, w^T Sigma w)) - f|^2
import numpy as np
from scipy.stats import norm
rng = np.random.default_rng(3)
n, N = 256, 40000
s2 = 2.0 / n
H = rng.standard_normal((N, n))
for _ in range(4):
    H = np.maximum(H @ (rng.standard_normal((n, n)) * np.sqrt(s2)), 0)
m = H.mean(0); e = m / np.linalg.norm(m); Q = np.eye(n) - np.outer(e, e)
Sig = np.cov(H.T, bias=True); SigP = Q @ Sig @ Q
Hp = H @ Q
gr = lambda mu, v: mu * norm.cdf(mu / np.sqrt(v)) + np.sqrt(v) * norm.pdf(mu / np.sqrt(v))   # E relu(N(mu,v))
bayes1, bayes2, gc, lev1 = [], [], [], []
for _ in range(20):
    w = rng.standard_normal(n) * np.sqrt(s2)
    a = H @ (np.outer(e, e) @ w)
    U = rng.standard_normal((6000, n)) * np.sqrt(s2) @ Q
    Z = a[:, None] + Hp @ U.T                         # (N, draws) pre-activations of the resampled rows
    F = np.maximum(Z, 0).mean(0)                       # exact neuron means for each resampled row
    q = np.einsum('di,ij,dj->d', U, SigP, U)
    # Bayes error given q: residual variance after flexible regression on q (degree-6 polynomial)
    Vq = np.vander((q - q.mean()) / q.std(), 7)
    resid = F - Vq @ np.linalg.lstsq(Vq, F, rcond=None)[0]
    bayes1.append(F.var()); bayes2.append(resid.var())
    # Gaussian closure using the exact mean/variance readouts of each resampled row
    mu = Z.mean(0); v = Z.var(0)
    gc.append(np.mean((gr(mu, v) - F)**2))
print(f"level-1 Bayes error  E Var(f | w.m)            = {np.mean(bayes1):.3e}")
print(f"level-2 Bayes error  E Var(f | w.m, w^T S w)    = {np.mean(bayes2):.3e}")
print(f"Gaussian-closure error (exact mean & variance) = {np.mean(gc):.3e}")
