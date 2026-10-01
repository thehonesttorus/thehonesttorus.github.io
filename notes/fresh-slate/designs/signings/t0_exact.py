"""T0: one copula step out of an exact Gaussian layer (layer 1) vs Monte Carlo cumulants of z_2."""
import sys, numpy as np, copula
n = int(sys.argv[1]) if len(sys.argv) > 1 else 12
rng = np.random.default_rng(1)
W = rng.standard_normal((2, n, n)) * np.sqrt(2 / n)
C1 = W[0].T @ W[0]; s = np.sqrt(np.diag(C1)); R = C1 / np.outer(s, s)
mu, st, dg = copula.layer_step(np.zeros(n), s, np.tile([1., 0, 0], (n, 1)), R, W[1])
# MC
N = int(float(sys.argv[2])) if len(sys.argv) > 2 else int(4e7)
S = np.zeros((4, n)); Sa = np.zeros(n); bs = 1 << 20; done = 0
while done < N:
    x = rng.standard_normal((bs, n))
    a = np.maximum(x @ W[0], 0); z = a @ W[1]
    Sa += a.sum(0)
    for k in range(4): S[k] += (z ** (k + 1)).sum(0)
    done += bs
m = S / done
mean = m[0]; v = m[1] - mean ** 2
c3 = m[2] - 3 * mean * m[1] + 2 * mean ** 3
c4c = m[3] - 4 * mean * m[2] + 6 * mean ** 2 * m[1] - 3 * mean ** 4
k4 = c4c - 3 * v ** 2
np.set_printoptions(precision=5, suppress=True, linewidth=150)
print("mean a1 err", np.abs(Sa / done - mu).max())
print("var  rel err", np.abs(dg['vz'] / v - 1).max())
print("k3 est", dg['k3'][:6]); print("k3 MC ", c3[:6])
print("k4 est", dg['k4'][:6]); print("k4 MC ", k4[:6])
for k, t in dg['terms'].items(): print(k, np.round(np.sqrt(np.mean(t ** 2)), 6))
print("k3 rms err", np.sqrt(np.mean((dg['k3'] - c3) ** 2)), "rms k3", np.sqrt(np.mean(c3 ** 2)))
print("k4 rms err", np.sqrt(np.mean((dg['k4'] - k4) ** 2)), "rms k4", np.sqrt(np.mean(k4 ** 2)))
