"""F0: exactness of the single-site formula where the Markov ansatz is exact (independent sites),
and the eps^2 law of the dropped two-site terms when sites are coupled."""
import sys, numpy as np
sys.path.insert(0, '.')
from mkv import mkv, relu_cov_offdiag
from scipy.integrate import quad
from scipy.stats import norm

def relu_cumulants(mu, s):
    m = [quad(lambda g: max(g, 0) ** k * norm.pdf(g, mu, s), mu - 14 * s, mu + 14 * s, points=[0], limit=200, epsabs=1e-14)[0] for k in range(5)]
    m1, m2, m3, m4 = m[1:]
    k2 = m2 - m1 ** 2; k3 = m3 - 3 * m2 * m1 + 2 * m1 ** 3
    k4 = m4 - 4 * m3 * m1 - 3 * m2 ** 2 + 12 * m2 * m1 ** 2 - 6 * m1 ** 4
    return k2, k3, k4

rng = np.random.default_rng(1); n = 8
W1 = np.diag(rng.uniform(0.5, 1.5, n)); W2 = rng.standard_normal((n, n)) * np.sqrt(2 / n)
W = np.stack([W1, W2, W2])
_, diag = mkv(W, w=4, return_all=True)
ex3 = np.zeros(n); ex4 = np.zeros(n)
for c in range(n):
    _, k3, k4 = relu_cumulants(0.0, W1[c, c])
    ex3 += W2[c] ** 3 * k3; ex4 += W2[c] ** 4 * k4
print("independent sites, layer-2 k3 rel err", np.abs(diag[1]['k3'] - ex3).max() / np.abs(ex3).max(),
      " k4 rel err", np.abs(diag[1]['k4'] - ex4).max() / np.abs(ex4).max())
# Mehler covariance vs MC
C = np.array([[1.0, 0.6], [0.6, 2.0]]); mu = np.array([0.3, -0.5])
x = rng.multivariate_normal(mu, C, 4_000_000); r = np.maximum(x, 0)
print("Mehler cov", relu_cov_offdiag(mu, C, 60)[0, 1], " MC", np.cov(r.T)[0, 1], "+-", r.std(0).prod() / 2000)
# coupled sites: layer-2 k3 error vs eps (truth by large MC)
for eps in [0.1, 0.2, 0.4]:
    A = np.eye(n) + eps * rng.standard_normal((n, n)) / np.sqrt(n)
    Wc = np.stack([A, W2, W2])
    _, dg = mkv(Wc, w=4, return_all=True)
    z = np.zeros((0, n)); tot = np.zeros(n); ks = []
    N = 0; S = np.zeros((4, n))
    for _ in range(20):
        xx = rng.standard_normal((1_000_000, n)); xx = np.concatenate([xx, -xx])
        z2 = np.maximum(xx @ A, 0) @ W2
        zc = z2
        S[0] += zc.sum(0); S[1] += (zc ** 2).sum(0); S[2] += (zc ** 3).sum(0); S[3] += (zc ** 4).sum(0); N += len(zc)
    m1, m2, m3, m4 = S / N
    k3 = m3 - 3 * m2 * m1 + 2 * m1 ** 3
    print(f"eps={eps}: k3 err rms {np.sqrt(np.mean((dg[1]['k3'] - k3) ** 2)):.2e}  (|k3| rms {np.sqrt(np.mean(k3**2)):.2e}; MC se ~{np.sqrt(15/N):.1e})")
