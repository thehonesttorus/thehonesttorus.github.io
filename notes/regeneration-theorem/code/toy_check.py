# Engine check for note XLIX on a toy with an exactly Gaussian layer: n = 6 sites, small cross-covariances.
# Compares the linked-cluster structure sum (all structures up to EMAX edges, no power-counting cut, Moebius-exact
# distinct-site sums) with a large Monte Carlo of the fourth cumulants of z' = W relu(z).
import sys, re, math, numpy as np
src = open(__file__.replace("toy_check.py", "regen4.py")).read()
head = src[src.index("def he("):src.index("def r2(")]
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
n = 6
EMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 5
N = int(float(sys.argv[3])) if len(sys.argv) > 3 else int(4e8)
SQ2PI = math.sqrt(2 * math.pi)
from scipy.special import ndtr
exec(head)
A = rng.normal(size=(n, n)) * 0.12
C = np.diag(rng.uniform(0.8, 1.6, n)) + 0.5 * (A + A.T)
np.fill_diagonal(C, np.diag(C))
assert np.all(np.linalg.eigvalsh(C) > 0)
mu = rng.normal(size=n) * 0.6
W = rng.normal(size=(3, n)) / math.sqrt(n)
Co = C.copy(); np.fill_diagonal(Co, 0.0)
F = site_factors(mu, np.diag(C).copy(), 0 * mu, 0 * mu, True)
mats = {"W": W, "Co": Co, "D21": 0 * Co, "K31": 0 * Co, "K22": 0 * Co}
res = {}
for name, rows in (("K31", "iiij"), ("K22", "iijj"), ("k4", "iiii")):
    tot = 0.0
    for st in structures(rows, hyper=False):
        Y = evaluate(rows, st, F, mats, outdiag=(name == "k4"))
        tot = tot + Y
    res[name] = tot
# Monte Carlo in chunks: central moments of y = W relu(z)
Lc = np.linalg.cholesky(C)
B = 2_000_000
m1 = np.zeros(3); S = {}
acc = {k: 0.0 for k in ("y", "yy", "y3y", "y2y2", "y4", "y2y", "y3")}
mom = None
sums = None
tot_n = 0
s1 = np.zeros(3); s2 = np.zeros((3, 3)); s3 = np.zeros((3, 3, 3)); s4 = np.zeros((3, 3, 3, 3))
for c in range(N // B):
    z = mu[None, :] + rng.normal(size=(B, n)) @ Lc.T
    y = np.maximum(z, 0.0) @ W.T
    s1 += y.sum(0); s2 += y.T @ y
    s3 += np.einsum("ni,nj,nk->ijk", y, y, y)
    s4 += np.einsum("ni,nj,nk,nl->ijkl", y, y, y, y)
    tot_n += B
m = s1 / tot_n; E2 = s2 / tot_n; E3 = s3 / tot_n; E4 = s4 / tot_n
# central moments
c2 = E2 - np.outer(m, m)
c3 = E3 - np.einsum("ij,k->ijk", E2, m) - np.einsum("ik,j->ijk", E2, m) - np.einsum("jk,i->ijk", E2, m) + 2 * np.einsum("i,j,k->ijk", m, m, m)
def c4f(i, j, k, l):
    mi, mj, mk, ml = m[i], m[j], m[k], m[l]
    return (E4[i, j, k, l] - mi * E3[j, k, l] - mj * E3[i, k, l] - mk * E3[i, j, l] - ml * E3[i, j, k]
            + mi * mj * E2[k, l] + mi * mk * E2[j, l] + mi * ml * E2[j, k] + mj * mk * E2[i, l] + mj * ml * E2[i, k]
            + mk * ml * E2[i, j] - 3 * mi * mj * mk * ml)
def k4(i, j, k, l):
    return c4f(i, j, k, l) - c2[i, j] * c2[k, l] - c2[i, k] * c2[j, l] - c2[i, l] * c2[j, k]
print(f"toy n={n}, EMAX={EMAX}, MC N={tot_n:.1e}")
for (i, j) in ((0, 1), (1, 2), (2, 0), (1, 0)):
    print(f"  K31[{i},{j}]: engine {res['K31'][i, j]:+.6f}  MC {k4(i, i, i, j):+.6f}   "
          f"K22[{i},{j}]: engine {res['K22'][i, j]:+.6f}  MC {k4(i, i, j, j):+.6f}")
for i in range(3):
    print(f"  k4[{i}]: engine {res['k4'][i]:+.6f}  MC {k4(i, i, i, i):+.6f}")
