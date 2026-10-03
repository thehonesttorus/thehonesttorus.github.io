# Neuron-side Cantor codes: row w_i -> (1{<w_i, x^(k)> > 0})_k for i.i.d. Gaussian inputs x^(k).
# (a) split depth of a pair ~ Geometric(theta/pi); E 2^{-k} = p/(1+p)
# (b) triple first split isolates i with prob (theta_ij + theta_ik - theta_jk)/(theta_ij+theta_ik+theta_jk)  (Gromov products)
# (c) number of level-k cylinders hit by the sphere = Cover's function C(k,d) = 2 sum_{i<d} binom(k-1,i)
# (d) collision mass sum_cells mu(cell)^2 = E (1 - Theta/pi)^k  ~ k^{-(d-1)}
import numpy as np
from math import comb
rng = np.random.default_rng(5)
d = 5
w = rng.standard_normal((3, d))
ang = lambda a, b: np.arccos(np.clip(a @ b / np.linalg.norm(a) / np.linalg.norm(b), -1, 1))
t01, t02, t12 = ang(w[0], w[1]), ang(w[0], w[2]), ang(w[1], w[2])
T, K = 200_000, 60
X = rng.standard_normal((T, K, d))
s = (X @ w.T > 0)                                       # (T, K, 3)
def first_diff(a, b):
    diff = s[:, :, a] != s[:, :, b]
    k = np.where(diff.any(1), diff.argmax(1) + 1, K + 1)
    return k
k01 = first_diff(0, 1); p = t01/np.pi
print(f"(a) pair: mean split depth {k01.mean():.4f} vs 1/p {1/p:.4f};  E 2^-k {np.mean(2.0**-k01):.4f} vs p/(1+p) {p/(1+p):.4f}")
# first split of the triple
iso = np.full(T, -1)
for t in range(K):
    b = s[:, t, :]; undec = iso < 0
    alone0 = (b[:, 0] != b[:, 1]) & (b[:, 1] == b[:, 2]); alone1 = (b[:, 1] != b[:, 0]) & (b[:, 0] == b[:, 2]); alone2 = (b[:, 2] != b[:, 0]) & (b[:, 0] == b[:, 1])
    iso[undec & alone0] = 0; iso[undec & alone1] = 1; iso[undec & alone2] = 2
tot = t01 + t02 + t12
pred = np.array([t01 + t02 - t12, t01 + t12 - t02, t02 + t12 - t01]) / tot
print(f"(b) triple first split: empirical {np.round(np.bincount(iso[iso>=0], minlength=3)/np.sum(iso>=0), 4)} vs Gromov {np.round(pred, 4)}")
# (c), (d): code the sphere by k random hyperplanes; collision mass averaged over 20 hyperplane draws
from math import gamma as G
M = 400_000
A, B = rng.standard_normal((1_000_000, d)), rng.standard_normal((1_000_000, d))
th = np.arccos(np.clip((A*B).sum(1)/np.linalg.norm(A, axis=1)/np.linalg.norm(B, axis=1), -1, 1))
Bd = np.sqrt(np.pi)*G((d-1)/2)/G(d/2)
for k in [4, 8, 16, 32, 64]:
    hits, coll = [], []
    for r in range(20):
        U = rng.standard_normal((M, d)); H = rng.standard_normal((k, d))
        codes = np.packbits(U @ H.T > 0, axis=1)
        codes = codes.view(np.dtype((np.void, codes.shape[1]))).ravel()
        _, cnt = np.unique(codes, return_counts=True)
        hits.append(len(cnt)); coll.append(np.sum((cnt/M)**2) - 1/M)     # unbiased for sum mu^2
    cover = 2*sum(comb(k-1, i) for i in range(d))
    asym = G(d-1)*np.pi**(d-1)/(Bd*k**(d-1))
    print(f"(c,d) k={k}: cylinders hit (max of 20) {max(hits)} <= Cover C(k,d)={cover}; collision mass {np.mean(coll):.5f} "
          f"vs E(1-Th/pi)^k {np.mean((1-th/np.pi)**k):.5f} (asymptotic {asym:.5f})")
