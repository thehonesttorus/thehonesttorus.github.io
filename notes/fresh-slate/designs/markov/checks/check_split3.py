"""Split kappa3(z_3) into fresh (layer-2 ReLU of a Gaussian z_2 with the true mean/cov) and old content;
compare each part with MKV's fresh-site and transported-site sums (n=64)."""
import sys, numpy as np
sys.path.insert(0, '.')
import mkv as M
from bake import weights
n = 64; N = int(2e7)
W = weights(n, 16, 1)[:3].astype(np.float64)
rng = np.random.default_rng(7)
def cum3(gen):
    S = np.zeros((3, n)); cnt = 0; S2 = np.zeros((n, n)); S1 = np.zeros(n)
    while cnt < N:
        z2, z3 = gen(200000)
        S[0] += z3.sum(0); S[1] += (z3**2).sum(0); S[2] += (z3**3).sum(0)
        S1 += z2.sum(0); S2 += z2.T @ z2; cnt += len(z3)
    m1, m2, m3 = S / cnt
    mu2 = S1 / cnt; C2 = S2 / cnt - np.outer(mu2, mu2)
    return m3 - 3*m2*m1 + 2*m1**3, mu2, C2
def real(b):
    x = rng.standard_normal((b, n)); z2 = np.maximum(x @ W[0], 0) @ W[1]
    return z2, np.maximum(z2, 0) @ W[2]
k3_real, mu2, C2 = cum3(real)
Lc = np.linalg.cholesky(C2)
def gauss(b):
    z2 = mu2 + rng.standard_normal((b, n)) @ Lc.T
    return z2, np.maximum(z2, 0) @ W[2]
k3_fresh, _, _ = cum3(gauss)
k3_old = k3_real - k3_fresh
# MKV parts at layer 3
mu = np.zeros(n); C = W[0].T @ W[0]
_, dg = M.mkv(W, w=4, var21='full', return_all=True)
d2 = dg[1]
m1, m2, p2 = M.relu_moments(d2['mu'], d2['v'], d2['k3'], d2['k4'])
coef2 = M.site_coeffs(d2['mu'], d2['v'])
# fresh: sites at layer 2 with MKV's C_2 replaced by the true C2 (to isolate the site formula)
fr3, _ = M.single_site(W[2], (C2 * p2[None, :]) @ W[2], coef2)
tot = dg[2]['k3']; old_mkv = tot - M.single_site(W[2], (np.diag(d2['v']) * 0 + 0) @ W[2] * 0 + 0, coef2)[0] * 0
r = lambda a, b: (np.sqrt(np.mean((a-b)**2))/np.sqrt(np.mean(b**2)), np.dot(a,b)/np.dot(b,b))
print("fresh: MKV-site(true C2) vs Gaussian-z2 truth: relerr %.3f slope %.3f" % r(fr3, k3_fresh))
print("rms fresh %.3e  rms old %.3e  rms total %.3e" % tuple(np.sqrt(np.mean(v**2)) for v in [k3_fresh, k3_old, k3_real]))
print("MKV total vs truth: relerr %.3f slope %.3f" % r(tot, k3_real))
print("MKV total - fresh(true C2) as old estimate vs truth old: relerr %.3f slope %.3f" % r(tot - fr3, k3_old))
