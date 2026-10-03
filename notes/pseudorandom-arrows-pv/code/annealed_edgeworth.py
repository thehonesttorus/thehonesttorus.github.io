# Does annealing the Edgeworth corrections (=> full-trace cumulant contractions) close the gap
# between Gaussian closure and the level-2 Bayes error?
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
Hp = H @ Q                                   # centred residual activations (Q m = 0)
xi = H @ e; xt = xi - xi.mean()
r2 = (Hp**2).sum(1)                          # |h'|^2
trS = np.trace(SigP)
# per-layer scalar contractions (computed once, O(N n)) -- the "full traces"
k3xi = np.mean(xt**3); k4xi = np.mean(xt**4) - 3 * np.mean(xt**2)**2
TrC3 = np.mean(xt * r2)                                            # sum_i k(xi,h_i,h_i)
cxh = (xt[:, None] * Hp).mean(0)
TrC4 = np.mean(xt**2 * r2) - np.mean(xt**2) * trS - 2 * cxh @ cxh  # sum_i k(xi,xi,h_i,h_i)
T4 = np.mean(r2**2) - trS**2 - 2 * np.sum(SigP**2)                 # sum_ik k4(h_i,h_i,h_k,h_k)
G = lambda mu, v: mu * norm.cdf(mu / np.sqrt(v)) + np.sqrt(v) * norm.pdf(mu / np.sqrt(v))
def d3(mu, v): return -(mu / v) * norm.pdf(mu / np.sqrt(v)) / np.sqrt(v)
def d4(mu, v): return ((mu**2 / v - 1) / v) * norm.pdf(mu / np.sqrt(v)) / np.sqrt(v)
def d6(mu, v):
    z = mu / np.sqrt(v); return (z**4 - 6 * z**2 + 3) / v**2 * norm.pdf(z) / np.sqrt(v)
gc, ed, ed_true, bayes = [], [], [], []
for _ in range(20):
    w = rng.standard_normal(n) * np.sqrt(s2)
    beta = w @ e
    U = rng.standard_normal((3000, n)) * np.sqrt(s2) @ Q
    Z = beta * xi[:, None] + Hp @ U.T
    F = np.maximum(Z, 0).mean(0)
    mu = Z.mean(0); v = Z.var(0)
    k3a = beta**3 * k3xi + 3 * beta * s2 * TrC3                     # annealed third cumulant
    k4a = beta**4 * k4xi + 6 * beta**2 * s2 * TrC4 + 3 * s2**2 * T4  # annealed fourth cumulant
    est_gc = G(mu, v)
    est_ed = est_gc + k3a / 6 * d3(mu, v) + k4a / 24 * d4(mu, v) + k3a**2 / 72 * d6(mu, v)
    Zc = Z - mu
    k3t = (Zc**3).mean(0); k4t = (Zc**4).mean(0) - 3 * v**2          # quenched (row-specific) cumulants
    est_et = est_gc + k3t / 6 * d3(mu, v) + k4t / 24 * d4(mu, v) + k3t**2 / 72 * d6(mu, v)
    q = np.einsum('di,ij,dj->d', U, SigP, U)
    Vq = np.vander((q - q.mean()) / q.std(), 7)
    bayes.append(np.var(F - Vq @ np.linalg.lstsq(Vq, F, rcond=None)[0]))
    gc.append(np.mean((est_gc - F)**2)); ed.append(np.mean((est_ed - F)**2)); ed_true.append(np.mean((est_et - F)**2))
print(f"Gaussian closure (exact mean, var)                   MSE = {np.mean(gc):.3e}")
print(f"+ annealed Edgeworth (full-trace cumulants, per layer) MSE = {np.mean(ed):.3e}")
print(f"+ quenched Edgeworth (exact per-row cumulants, oracle)  MSE = {np.mean(ed_true):.3e}")
print(f"level-2 Bayes error (reference)                         = {np.mean(bayes):.3e}")
