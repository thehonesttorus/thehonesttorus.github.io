# Independent numerical check of the null-source identities (THEORY.md eq. (2) and the degree-p generalisation).
# For a degree-p homogeneous F and Z ~ N(mu, Sigma) in 2D:  L_k[T]F = (1/k!) E[F(Z) T : H_k(Z)], with the Gaussian Hermite
# tensors H_k = (-1)^k nabla^k phi / phi:  y = Q (z - mu), Q = Sigma^-1,
#   H2 = y y - Q,   H3 = y y y - 3 sym(y Q),   H4 = y y y y - 6 sym(y y Q) + 3 sym(Q Q).
# Claim (degree p):  L4[Sigma.G] + (1/4) L3[mu.G] - ((p - 2)/12) L2[G] = 0:  null tuple (G/12, mu.G/4, Sigma.G) for p = 1,
# (0, mu.G/4, Sigma.G) for p = 2. Integrals by sector-wise polar quadrature (F smooth inside sectors cut by gate rays).
import itertools, math, numpy as np
from scipy import integrate
mu = np.array([0.3, -0.2]); Sig = np.array([[1.0, 0.4], [0.4, 0.7]]); Gm = np.array([[0.5, -0.3], [-0.3, 0.8]])
Q = np.linalg.inv(Sig); norm = 1.0 / (2 * np.pi * math.sqrt(np.linalg.det(Sig)))
perm4 = list(itertools.permutations(range(4)))
S4 = np.zeros((2,) * 4); M3 = np.zeros((2,) * 3)
for idx in itertools.product(range(2), repeat=4):
    S4[idx] = sum(Sig[idx[p[0]], idx[p[1]]] * Gm[idx[p[2]], idx[p[3]]] for p in perm4) / 24
for i, j, k in itertools.product(range(2), repeat=3):
    M3[i, j, k] = (mu[i] * Gm[j, k] + mu[j] * Gm[i, k] + mu[k] * Gm[i, j]) / 3
def sym(T):
    k = T.ndim; return sum(np.transpose(T, p) for p in itertools.permutations(range(k))) / math.factorial(k)
QQs = sym(np.einsum("ij,kl->ijkl", Q, Q))
def herm_pair(x, y_):
    z = np.array([x, y_]); yv = Q @ (z - mu); ph = norm * math.exp(-0.5 * (z - mu) @ Q @ (z - mu))
    yy = np.outer(yv, yv)
    H2 = yy - Q
    H3 = np.einsum("i,j,k->ijk", yv, yv, yv) - 3 * sym(np.einsum("i,jk->ijk", yv, Q))
    H4 = np.einsum("i,j,k,l->ijkl", yv, yv, yv, yv) - 6 * sym(np.einsum("i,j,kl->ijkl", yv, yv, Q)) + 3 * QQs
    return ph * np.sum(Gm * H2) / 2, ph * np.sum(M3 * H3) / 6, ph * np.sum(S4 * H4) / 24
relu = lambda v: np.maximum(v, 0.0)
W1 = np.array([[1.0, -0.6], [0.4, 0.9], [-0.8, 0.5]]); W2 = np.array([[0.5, -0.4, 0.8], [-0.3, 0.9, 0.2]]); w2 = np.array([0.7, -1.1, 0.9])
def F1(x, y_):   # degree 1: two-layer bias-free ReLU net with a signed readout
    h = relu(W1 @ np.array([x, y_])); return float(w2 @ h - 0.3 * relu(W2 @ h).sum())
def F2(x, y_):   # degree 2: products of hidden units
    h = relu(W1 @ np.array([x, y_])); g = relu(W2 @ h); return float(g[0] * g[1] + h[0] * h[2])
def rays(F):     # kink angles: where the second difference of F on the unit circle jumps
    th = np.linspace(0, 2 * np.pi, 400001)[:-1]; v = np.array([F(math.cos(t), math.sin(t)) for t in th])
    d2v = np.abs(np.diff(np.r_[v[-1], v, v[0]], 2)); thr = 1e3 * np.median(d2v) + 1e-14
    ks = []
    for i in np.where(d2v > thr)[0]:
        t = th[i]
        if not ks or t - ks[-1] > 1e-3: ks.append(t)
    return [0.0] + ks + [2 * np.pi]
for name, F, p in (("degree 1", F1, 1), ("degree 2", F2, 2)):
    ks = rays(F); acc = np.zeros(3)
    for a, b in zip(ks[:-1], ks[1:]):
        if b - a < 1e-12: continue
        for c in range(3):
            val, _ = integrate.dblquad(lambda r, t: F(r * math.cos(t), r * math.sin(t)) * herm_pair(r * math.cos(t), r * math.sin(t))[c] * r,
                                       a, b, 0.0, 14.0, epsabs=1e-13, epsrel=1e-11)
            acc[c] += val
    L2, L3, L4 = acc
    nul = L4 + L3 / 4 - (p - 2) / 12 * L2
    other = L4 + L3 / 4 + L2 / 12 if p == 2 else L4 + L3 / 4     # the other degree's tuple, as a control
    print(f"{name}: L2[G] {L2:+.10f}  L3[mu.G] {L3:+.10f}  L4[Sigma.G] {L4:+.10f} | null combination {nul:+.2e} | "
          f"control (other degree's tuple) {other:+.2e}  [{len(ks) - 2} kink rays]", flush=True)
