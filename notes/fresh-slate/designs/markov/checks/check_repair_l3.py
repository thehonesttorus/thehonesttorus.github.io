"""Test of the MKV-2 repair at layer 3: a layer-1 site's emission crosses the layer-2 ReLU as a function of
its own latent xi: G(xi) = sum_b W3_bj [psi'_b P_sb X + psi''_b (c_sb xi P_sb X + 1/2 P_sb^2 (X^2 - E X^2))],
He1 part folded into the link c.  Compare layer-3 kappa3 slope vs MC truth: linear transport vs repaired."""
import sys, numpy as np
sys.path.insert(0, '.')
import mkv as M
from bake import weights
n = int(sys.argv[1]); N = int(float(sys.argv[2]))
W = weights(n, 16, 1)[:3].astype(np.float64); rng = np.random.default_rng(9)
S = np.zeros((3, n)); cnt = 0
while cnt < N:
    x = rng.standard_normal((100000, n)); x = np.concatenate([x, -x])
    z3 = np.maximum(np.maximum(x @ W[0], 0) @ W[1], 0) @ W[2]
    S[0] += z3.sum(0); S[1] += (z3**2).sum(0); S[2] += (z3**3).sum(0); cnt += len(z3)
m1, m2, m3 = S / cnt; k3t = m3 - 3*m2*m1 + 2*m1**3

_, dg = M.mkv(W, w=4, var21='full', return_all=True)
# layer-1 sites
C1 = W[0].T @ W[0]; s1 = np.sqrt(np.diag(C1)); mu1 = np.zeros(n)
ma, _, p1 = M.relu_moments(mu1, s1**2)
coef1 = M.site_coeffs(mu1, s1**2)
# layer-2 state from MKV
d2 = dg[1]; mu2, v2 = d2['mu'], d2['v']; sd2 = np.sqrt(v2)
_, _, p2 = M.relu_moments(mu2, v2, d2['k3'], d2['k4'])
dens2 = M.pdf(mu2 / sd2) / sd2
P12 = W[1]; K12 = (C1 * p1[None, :]) @ W[1]          # site s -> z_2b ; Cov(g_s, z_2b)
c12 = K12 / s1[:, None]                               # standardised link
# site functions on quadrature nodes (standardised xi), per site s
xg, wg = np.polynomial.hermite_e.hermegauss(200); wg = wg / wg.sum()
# kink-robust: use fine Gauss-Legendre pieces instead
lo, hi = -12.0, 12.0
u0 = np.zeros(n)  # mu1 = 0 -> kink at 0
GLx, GLw = np.polynomial.legendre.leggauss(150)
nodes = np.concatenate([(lo + 0) / 2 + (0 - lo) / 2 * GLx, (0 + hi) / 2 + (hi - 0) / 2 * GLx])
wts = np.concatenate([(0 - lo) / 2 * GLw, (hi - 0) / 2 * GLw]) * M.pdf(nodes)
xi = nodes[None, :]
X = np.maximum(s1[:, None] * xi, 0) - ma[:, None] - p1[:, None] * s1[:, None] * xi
EX2 = (wts * X**2).sum(1)
U = np.stack([X, xi * X, X**2 - EX2[:, None]])        # (3, n_sites, Q)
lam = (wts * U * xi).sum(-1)                          # He1 coefficients E[u_i xi]  (3, n)
Ut = U - lam[:, :, None] * xi[None]                   # project out He1
M3 = np.einsum('q,isq,jsq,ksq->sijk', wts, Ut, Ut, Ut)
M2x = np.einsum('q,isq,jsq->sij', wts, Ut * xi[None], Ut)
Mh2 = np.einsum('q,isq->si', wts, Ut * (xi**2 - 1)[None])
# channel coefficients into z_3j
Pp = (P12 * p2[None, :]) @ W[2]
R = (P12 * c12 * dens2[None, :]) @ W[2]
Q = 0.5 * (P12**2 * dens2[None, :]) @ W[2]
K13 = (K12 * p2[None, :]) @ W[2] / s1[:, None]       # standardised linear link to z_3j
def old_k3(a, c):
    t3 = np.einsum('isj,jsk... ->', np.zeros((1,1,1)), np.zeros((1,1,1))) if False else None
    A = a  # (3, s, j)
    e3 = np.einsum('sijk,isj,jsj2->sj2', M3, A, A) if False else None
    E3 = np.einsum('sijk,isn,jsn,ksn->sn', M3, A, A, A)
    E2 = np.einsum('sij,isn,jsn->sn', M2x, A, A)
    E1 = np.einsum('si,isn->sn', Mh2, A)
    return (E3 + 3 * c * E2 + 3 * c**2 * E1).sum(0)
lin = old_k3(np.stack([Pp, 0 * R, 0 * Q]), K13)
rep_c = K13 + lam[1][:, None] * R + lam[2][:, None] * Q
rep = old_k3(np.stack([Pp, R, Q]), rep_c)
# fresh layer-2 sites (MKV) = MKV total minus MKV's linear old part
fresh = dg[2]['k3'] - lin
for name, est in [('linear transport', fresh + lin), ('repaired', fresh + rep)]:
    print(f"n={n} {name:17s} slope {np.dot(est,k3t)/np.dot(k3t,k3t):.3f} relerr {np.sqrt(np.mean((est-k3t)**2)/np.mean(k3t**2)):.3f}", flush=True)
