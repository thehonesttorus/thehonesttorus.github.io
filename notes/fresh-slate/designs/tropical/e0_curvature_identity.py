"""E0: exactness check of the tropical-curvature identity (I2) at toy scale.

E[a_{l,j}] = sum_{k<=l} sum_m E[ delta(z_{k,m}) |grad_x z_{k,m}|^2 S_{(l,j)<-(k,m)} ]

Both sides are integrated along random lines x = x_perp + t v (t ~ N(0,1), v a random unit vector):
LHS: the 1-D piecewise-linear restriction integrated against phi(t) (trapezoid on a fine grid);
RHS: E[delta(z) H] = E_{x_perp} sum_{roots t*} phi(t*) H(t*) / |v . grad z(t*)|  (roots by linear
interpolation inside the grid cell: exact unless a second kink falls in the same cell).
The lines are shared by both sides, so the comparison is paired; the per-line identity does NOT hold
(only its average), so the residual is MC noise over x_perp, reported with its standard error.
Usage: python e0_curvature_identity.py seed n L M
"""
import numpy as np, sys
seed, n, L, M = [int(a) for a in (sys.argv[1:] + ['0', '8', '3', '20000'][len(sys.argv) - 1:])][:4]
rng = np.random.default_rng(seed)
W = [rng.standard_normal((n, n)) * np.sqrt(2 / n) for _ in range(L)]
T = np.linspace(-9, 9, 3601); h = T[1] - T[0]
phi = lambda t: np.exp(-t * t / 2) / np.sqrt(2 * np.pi)
wq = phi(T) * h; wq[0] *= .5; wq[-1] *= .5

def forward(X):
    zs = []; a = X
    for l in range(L):
        z = a @ W[l]; zs.append(z); a = np.maximum(z, 0)
    return zs

def batch_jac(X):
    """X (B,n): per layer (z (B,n), grad G (B,n,n) [input i, neuron m], gate (B,n))"""
    B = X.shape[0]; a = X; J = np.broadcast_to(np.eye(n), (B, n, n)).copy(); out = []
    for l in range(L):
        z = a @ W[l]; G = J @ W[l]; g = (z > 0).astype(float)
        out.append((z, G, g)); a = np.maximum(z, 0); J = G * g[:, None, :]
    return out

lhs_s = np.zeros((L, n)); rhs_s = np.zeros((L, n)); d2 = np.zeros((L, n))
CH = 200
for c0 in range(0, M, CH):
    V = rng.standard_normal((CH, n)); V /= np.linalg.norm(V, axis=1, keepdims=True)
    X0 = rng.standard_normal((CH, n)); X0 -= (X0 * V).sum(1, keepdims=True) * V
    Xg = X0[:, None, :] + T[None, :, None] * V[:, None, :]          # (CH, nt, n)
    zs = forward(Xg.reshape(-1, n))
    lhs_line = np.stack([np.einsum('t,ctn->cn', wq, np.maximum(z.reshape(CH, len(T), n), 0)) for z in zs], 1)  # (CH,L,n)
    rhs_line = np.zeros((CH, L, n))
    for k in range(L):
        zk = zs[k].reshape(CH, len(T), n)
        ci, ti, mi = np.nonzero(np.sign(zk[:, :-1, :]) * np.sign(zk[:, 1:, :]) < 0)
        z0 = zk[ci, ti, mi]; z1 = zk[ci, ti + 1, mi]
        tr = T[ti] - z0 * h / (z1 - z0)
        Xr = X0[ci] + tr[:, None] * V[ci]
        st = batch_jac(Xr)
        grad = st[k][1][np.arange(len(ci)), :, mi]                       # (R, n)
        base = np.nan_to_num(phi(tr) * (grad * grad).sum(1) / np.abs((V[ci] * grad).sum(1)))
        np.add.at(rhs_line, (ci, k, mi), base)
        sens = np.zeros((len(ci), n)); sens[np.arange(len(ci)), mi] = 1.0
        for l in range(k + 1, L):
            sens = (sens @ W[l]) * st[l][2]
            np.add.at(rhs_line, (ci, l), base[:, None] * sens)
    lhs_s += lhs_line.sum(0); rhs_s += rhs_line.sum(0); d2 += ((lhs_line - rhs_line) ** 2).sum(0)
lhs = lhs_s / M; rhs = rhs_s / M; se = np.sqrt(d2 / M - (lhs - rhs) ** 2) / np.sqrt(M)
exact1 = np.linalg.norm(W[0], axis=0) / np.sqrt(2 * np.pi)
print(f"n={n} L={L} M={M}")
print("layer-1 closed form |w|/sqrt(2pi) vs LHS: max dev %.2e ; vs RHS: max dev %.2e" % (np.abs(lhs[0] - exact1).max(), np.abs(rhs[0] - exact1).max()))
for l in range(L):
    print("layer", l + 1, "max|LHS-RHS| %.2e  max z-score %.2f  (max|LHS| %.3f)" % (np.abs(lhs[l] - rhs[l]).max(), (np.abs(lhs[l] - rhs[l]) / se[l]).max(), np.abs(lhs[l]).max()))
