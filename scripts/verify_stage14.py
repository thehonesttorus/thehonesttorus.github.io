"""Checks for the stage-14 note (lattices, leaves and ellipsoids).

  1. Klartag's terminal states are Voronoi-perfect; local optima are eutactic (John). A2 against Z^2.
  2. The kink chaos law: the k-th chaos energy of relu(r + xi) is He_{k-2}(r)^2 phi(r)^2 / k!; it sums to
     Var relu - Phi(r)^2; for large |r| it sits at k ~ r^2; at r = 0 its tail beyond K decays like K^{-3/2}.
  3. The maximal resolved probe about a fixed centre is the Loewner ellipsoid of the dual wall points; its whitened contact
     normals form a tight frame (John's decomposition), sum of weights = n.
  4. Homogeneous walls: the dilation (dm, dA) = (m, 2A) is tangent to every contact constraint; affine walls break it.
  5. Conic Steiner formula in a network cell C = W^-1 R^n_+ (W Gaussian): with v_k the law of the face dimension of the
     projection, E|P_C g|^2 = sum k v_k, E|P_C g|^4 = sum v_k k(k+2), E|P_C g|^2 |P_C° g|^2 = sum v_k k(n-k).
  6. The sticky (Klartag) evolution with walls as obstacles: Ito accounting of log det A along a simulated path.
"""
import numpy as np
from math import pi, sqrt, log, lgamma, exp
from scipy.stats import norm
from scipy.optimize import nnls

phi, Phi = norm.pdf, norm.cdf
rng = np.random.default_rng(14)

print("== 1. Perfect and eutactic: A2 against Z^2")
for name, B in [("A2", np.array([[1.0, 0.0], [0.5, sqrt(3) / 2]])), ("Z2", np.eye(2))]:
    pts = np.array([i * B[0] + j * B[1] for i in range(-3, 4) for j in range(-3, 4) if (i, j) != (0, 0)])
    nr = np.einsum("ij,ij->i", pts, pts)
    mv = pts[np.isclose(nr, nr.min())]
    sym = np.array([[x[0] ** 2, x[0] * x[1], x[1] ** 2] for x in mv])
    S = sum(np.outer(x, x) for x in mv)
    print("  %s: %d minimal vectors; rank of {x x^T} = %d (perfect iff 3); sum x x^T = %s (eutactic iff prop. to I)"
          % (name, len(mv), np.linalg.matrix_rank(sym), np.round(S.ravel(), 6)))

print("== 2. The kink chaos law")


def kink_energies(r, K):
    h = np.zeros(K + 1)
    h[0], h[1] = 1.0, r
    for j in range(1, K):
        h[j + 1] = (r * h[j] - sqrt(j) * h[j - 1]) / sqrt(j + 1)
    k = np.arange(2, K + 3)
    return k, h[: K + 1] ** 2 * phi(r) ** 2 / (k * (k - 1))


m1 = lambda r: phi(r) + r * Phi(r)
m2 = lambda r: (1 + r * r) * Phi(r) + r * phi(r)
for r in (0.0, 1.0, 3.0, 5.0, 8.0):
    k, E = kink_energies(r, 20000)
    tot = E.sum()
    print("  r=%3.1f: sum_k E_k = %.6e ; Var relu - Phi^2 = %.6e ; mean chaos order %.2f (r^2 = %.1f), median %d"
          % (r, tot, m2(r) - m1(r) ** 2 - Phi(r) ** 2, (k * E).sum() / tot, r * r, k[np.searchsorted(np.cumsum(E), tot / 2)]))
k, E = kink_energies(0.0, 200000)
tails = [(K, E[k > K].sum()) for K in (100, 1000, 10000)]
print("  r=0 tail beyond K: " + ", ".join("K=%d %.3e" % t for t in tails) +
      " ; slope %.3f, %.3f (K^-3/2 -> -1.5)" % (log(tails[1][1] / tails[0][1]) / log(10), log(tails[2][1] / tails[1][1]) / log(10)))
for r in (4.0, 8.0, 12.0):
    k, E = kink_energies(r, 40000)
    tp = r * r / 4
    sh = lambda K: E[k <= K].sum() / E.sum()
    tail = lambda K: E[k > K].sum()
    print("  r=%.0f: share below the turning point k = r^2/4: %.3f ; share in chaos <= r^2/8: %.4f ; tail exponent beyond 4 r^2: %.3f"
          % (r, sh(tp), sh(tp / 2), log(tail(16 * tp) / tail(8 * tp)) / log(2)))
for r in (3.0, 5.0, 8.0):
    k, E = kink_energies(r, 20000)
    print("  r=%.0f: share of the kink energy in chaos <= K for K = r^2/2, r^2, 2r^2: %.3f, %.3f, %.3f"
          % (r, E[k <= 0.5 * r * r].sum() / E.sum(), E[k <= r * r].sum() / E.sum(), E[k <= 2 * r * r].sum() / E.sum()))

print("== 3. Maximal resolved probe = Loewner ellipsoid of the dual wall points; tight frame")
n = 4
Wd = rng.normal(size=(9, n))
b = rng.normal(size=9) * 0.5
m = rng.normal(size=n) * 0.3
sgn = np.sign(Wd @ m - b)
d = np.abs(Wd @ m - b)
rho = 1.0
X = rho * Wd / d[:, None]
P = np.vstack([X, -X])
u = np.ones(len(P)) / len(P)
for it in range(20000):
    Minv = np.linalg.inv(P.T @ (u[:, None] * P))
    M = np.einsum("ij,jk,ik->i", P, Minv, P)
    j = M.argmax()
    step = (M[j] - n) / (n * (M[j] - 1))
    u = (1 - step) * u
    u[j] += step
    if M[j] - n < 1e-10:
        break
A = np.linalg.inv(P.T @ (u[:, None] * P)) / n
cons = np.einsum("ij,jk,ik->i", X, A, X)
act = cons > 1 - 1e-6
lam = n * (u[: len(X)] + u[len(X):])
U = X @ np.linalg.cholesky(A)
print("  active walls %d of %d ; max_a x_a^T A x_a = %.6f (resolved iff <= 1)" % (act.sum(), len(X), cons.max()))
print("  || sum_a lam_a u_a u_a^T - I || = %.2e ; sum lam = %.4f (n = %d)"
      % (np.linalg.norm(U.T @ (lam[:, None] * U) - np.eye(n)), lam.sum(), n))
dist = np.abs(Wd @ m - b) / np.sqrt(np.einsum("ij,jk,ik->i", Wd, A, Wd))
print("  Mahalanobis distance of each wall from the probe centre (contacts = rho = 1):", np.round(np.sort(dist)[:6], 4))

print("== 4. Homogeneous walls keep the dilation free")
Wh = rng.normal(size=(6, 3))
mh = rng.normal(size=3)
Ah_ = np.eye(3) * 0.2
for walls, bb in (("homogeneous", np.zeros(6)), ("affine", rng.normal(size=6))):
    g = (Wh @ mh - bb) ** 2 - np.einsum("ij,jk,ik->i", Wh, Ah_, Wh)
    dm, dA = mh, 2 * Ah_
    dg = 2 * (Wh @ mh - bb) * (Wh @ dm) - np.einsum("ij,jk,ik->i", Wh, dA, Wh)
    print("  %-11s walls: d g_a along the dilation minus 2 g_a = %s" % (walls, np.round(dg - 2 * g, 6)))

print("== 5. Conic Steiner formula in a layer-1 cell C = W^-1 R^n_+")
n, N = 10, 20000
W = rng.normal(size=(n, n)) * sqrt(2 / n)
Cgen = np.linalg.inv(W)
G = rng.normal(size=(N, n))
a2, b2, kdim = np.empty(N), np.empty(N), np.empty(N, dtype=int)
for i in range(N):
    y, _ = nnls(Cgen, G[i])
    p = Cgen @ y
    a2[i] = p @ p
    b2[i] = (G[i] - p) @ (G[i] - p)
    kdim[i] = (y > 1e-12).sum()
v = np.bincount(kdim, minlength=n + 1) / N
kk = np.arange(n + 1)
se = lambda z: z.std() / sqrt(N)
print("  face-dimension law v_k:", np.round(v, 4), " statistical dimension sum k v_k = %.4f (orthogonal W: %.1f)" % ((kk * v).sum(), n / 2))
print("  E|P g|^2   = %.4f (se %.4f) ; sum k v_k        = %.4f" % (a2.mean(), se(a2), (kk * v).sum()))
print("  E|P g|^4   = %.3f (se %.3f) ; sum v_k k(k+2)   = %.3f" % ((a2 ** 2).mean(), se(a2 ** 2), (v * kk * (kk + 2)).sum()))
print("  E|P g|^2|P° g|^2 = %.3f (se %.3f) ; sum v_k k(n-k) = %.3f" % ((a2 * b2).mean(), se(a2 * b2), (v * kk * (n - kk)).sum()))

print("== 6. Sticky (Klartag) evolution with walls as obstacles: Ito accounting of log det A (A = probe precision)")
n = 3
nw = 14
Ww = rng.normal(size=(nw, n))
bw = rng.normal(size=nw) * 0.6
mc = np.zeros(n)
cw = (Ww @ mc - bw) ** 2  # probe {(x-m)^T A (x-m) < 1} avoids wall a iff w^T A^-1 w <= c_a
a0 = 2 * np.max(np.einsum("ij,ij->i", Ww, Ww) / cw)
iu = [(i, j) for i in range(n) for j in range(i, n)]


def to_mat(v):
    M = np.zeros((n, n))
    for k, (i, j) in enumerate(iu):
        if i == j:
            M[i, i] = v[k]
        else:
            M[i, j] = M[j, i] = v[k] / sqrt(2)
    return M


def to_vec(M):
    return np.array([M[i, j] if i == j else sqrt(2) * M[i, j] for (i, j) in iu])


sig, T, dt, paths = 0.25 * a0, 2.0, 2e-3, 200
final_ld, pred, ncont = [], [], []
for p_ in range(paths):
    A = a0 * np.eye(n)
    Y = []
    drift = 0.0
    for step in range(int(T / dt)):
        if Y:
            Q, _ = np.linalg.qr(np.array([to_vec(np.outer(y, y)) for y in Y]).T)
            if Q.shape[1] >= len(iu):
                break
            Pi = np.eye(len(iu)) - Q @ Q.T
        else:
            Pi = np.eye(len(iu))
        lam_, UA = np.linalg.eigh(A)
        delta = 0.0
        for i in range(n):
            for j in range(n):
                c = Pi @ to_vec((np.outer(UA[:, i], UA[:, j]) + np.outer(UA[:, j], UA[:, i])) / 2)
                delta += (c @ c) / (lam_[i] * lam_[j])
        drift += -0.5 * sig ** 2 * delta * dt
        A = A + sig * to_mat(Pi @ rng.normal(size=len(iu))) * sqrt(dt)
        Ainv = np.linalg.inv(A)
        g = np.einsum("ij,jk,ik->i", Ww, Ainv, Ww) - cw
        for a in np.where(g >= 0)[0]:
            y = Ainv @ Ww[a]
            if all(abs(np.dot(y, yy) / (np.linalg.norm(y) * np.linalg.norm(yy))) < 1 - 1e-9 for yy in Y):
                Y.append(y)
    final_ld.append(np.linalg.slogdet(A)[1])
    pred.append(n * log(a0) + drift)
    ncont.append(len(Y))
final_ld, pred = np.array(final_ld), np.array(pred)
print("  E log det A_T = %.4f (se %.4f) ; log det A_0 - (1/2) E int delta = %.4f (se %.4f) ; log det A_0 = %.4f"
      % (final_ld.mean(), final_ld.std() / sqrt(paths), pred.mean(), pred.std() / sqrt(paths), n * log(a0)))
print("  contacts at T: mean %.2f, max %d ; dim Sym_3 = 6 (the process freezes when the contact y y^T span it)"
      % (np.mean(ncont), max(ncont)))
