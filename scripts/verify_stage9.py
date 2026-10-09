"""Numerical checks of the stage-9 identities.
  python scripts/verify_stage9.py
1. Index shift d_k' = d_{k+1} (Lemma).  2. Price/Stein: dM/dm_a = E[theta_a relu_b], dM/dS_ab = E[theta_a theta_b],
d^2 M/dm_a dm_b = dM/dS_ab (one-layer heat equation).  3. Theorem (defect = curvature x Gram) on a tiny network:
delta(0) = e - y.grad e - Lap e by finite differences in y, against -<Hess_(mu,M) G, Gamma> with G the chain from
layer 2 as a function of the first-layer moments (Hessian by finite differences) and Gamma from the closed forms."""
import sys, os, math, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from whest.relu_gauss import hermite_relu, relu_moments, phi
from scipy.special import ndtr
rng = np.random.default_rng(1)
# 1. index shift
al = np.linspace(-3, 3, 25); K = 10; eps = 1e-5
d = hermite_relu(al, K + 1); dp = (hermite_relu(al + eps, K + 1) - hermite_relu(al - eps, K + 1)) / (2 * eps)
print(f"1. index shift: max|d_k' - d_(k+1)| (k<={K}) = {np.max(np.abs(dp[:K+1] - d[1:K+2])):.1e}")
# 2. Price / Stein on a bivariate Gaussian
m = np.array([0.3, -0.5]); S = np.array([[1.3, 0.4], [0.4, 0.8]])
def Mfun(m, S):
    mu, Kh, sig, alp = relu_moments(m, S, 40); return Kh + np.outer(mu, mu)
def fd(f, x, i, j=None, h=1e-4):
    e = np.zeros_like(x); e[i] = h
    if j is None: return (f(x + e) - f(x - e)) / (2 * h)
    e2 = np.zeros_like(x); e2[j] = h
    return (f(x + e + e2) - f(x + e - e2) - f(x - e + e2) + f(x - e - e2)) / (4 * h * h)
N = 4_000_000; Lc = np.linalg.cholesky(S); z = m[:, None] + Lc @ rng.standard_normal((2, N)); r = np.maximum(z, 0); th = (z > 0)
dM_dma = fd(lambda mm: Mfun(mm, S)[0, 1], m, 0)
dM_dSab = (Mfun(m, S + 1e-4 * np.array([[0, 1], [1, 0]]))[0, 1] - Mfun(m, S - 1e-4 * np.array([[0, 1], [1, 0]]))[0, 1]) / 2e-4
d2M = fd(lambda mm: Mfun(mm, S)[0, 1], m, 0, 1)
print(f"2. dM/dm_a = {dM_dma:.5f}  MC E[theta_a relu_b] = {np.mean(th[0] * r[1]):.5f} | dM/dS_ab = {dM_dSab:.5f}  MC E[theta_a theta_b] = {np.mean(th[0] * th[1]):.5f}  d2M/dm_a dm_b = {d2M:.5f}")
# 3. Gram theorem on a tiny network (Gaussian closure), n small so the Hessian of G in (mu, M) is cheap
n, L = 6, 3; W = rng.standard_normal((L, n, n)) * np.sqrt(2.0 / n)
def layer(mu, C, l):
    mz = W[l] @ mu; Cz = W[l] @ C @ W[l].T; m2, K2, _, _ = relu_moments(mz, Cz, 60); return m2, K2
def chain(y):
    mu, C = y.copy(), np.eye(n)
    for l in range(L): mu, C = layer(mu, C, l)
    return mu
def theta(y):                      # exact first-layer moments (mu_1, M_1) for input N(y, I)
    mu, Kh = layer(y, np.eye(n), 0); return mu, Kh + np.outer(mu, mu)
def G(mu, M):                      # chain from layer 2 as a function of (mu_1, M_1)
    C = M - np.outer(mu, mu)
    for l in range(1, L): mu, C = layer(mu, C, l)
    return mu
h = 1e-3; y0 = np.zeros(n); e0 = chain(y0)
grad = np.array([fd(chain, y0, i, h=h) for i in range(n)])            # (n, n_out)
lap = sum((chain(y0 + h * np.eye(n)[i]) + chain(y0 - h * np.eye(n)[i]) - 2 * e0) / h ** 2 for i in range(n))
delta_fd = e0 - lap                                                   # y.grad e = 0 at y = 0
# Gram from closed forms at y=0
mu1, M1 = theta(y0); mz = W[0] @ y0; Sz = W[0] @ W[0].T; sig = np.sqrt(np.diag(Sz)); alp = mz / sig
g = ndtr(alp); dK = hermite_relu(alp, 60); rho = Sz / np.outer(sig, sig); np.fill_diagonal(rho, 0)
T = np.zeros((n, n))                                                  # T_ab = dM_ab/dm_a = sigma_b sum rho^k d_{k+1}(a) d_k(b)/k!
for k in range(59): T += np.outer(dK[k + 1], dK[k]) * rho ** k / math.factorial(k)
T *= sig[None, :]
# diagonal: d M_aa / d m_a = 2 E[theta_a relu_a] = 2 sigma_a (phi(alpha) + alpha Phi(alpha))... use finite difference of the exact diagonal
def Mdiag(mm, a):
    s = sig[a]; al_ = mm[a] / s; return s * s * ((1 + al_ * al_) * ndtr(al_) + al_ * phi(al_))
for a in range(n):
    T[a, a] = (Mdiag(mz + 1e-5 * np.eye(n)[a], a) - Mdiag(mz - 1e-5 * np.eye(n)[a], a)) / 2e-5
# moment index set: A = (mu_a) for a<n, then (M_ab) for a<=b ; gradients in y: grad mu_a = g_a w_a ; grad M_ab = T_ab w_a + T_ba w_b (a!=b), grad M_aa = T_aa w_a
idx = [("mu", a, a) for a in range(n)] + [("M", a, b) for a in range(n) for b in range(a, n)]
grads = []
for kind, a, b in idx:
    if kind == "mu": grads.append(g[a] * W[0][a])
    elif a == b: grads.append(T[a, a] * W[0][a])
    else: grads.append(T[a, b] * W[0][a] + T[b, a] * W[0][b])
grads = np.array(grads); Gamma = grads @ grads.T
# check the gradient closed forms against finite differences of theta
grad_fd = []
for kind, a, b in idx:
    f = (lambda yy, a=a: theta(yy)[0][a]) if kind == "mu" else (lambda yy, a=a, b=b: theta(yy)[1][a, b])
    grad_fd.append(np.array([fd(f, y0, i, h=h) for i in range(n)]))
print(f"3. first-layer moment gradients: closed form vs finite differences, max abs diff {np.max(np.abs(np.array(grad_fd) - grads)):.1e} (scale {np.max(np.abs(grads)):.2f})")
# Hessian of G in the moment coordinates (symmetric M perturbed symmetrically)
def Gvec(v):
    mu = v[:n]; M = M1.copy(); c = n
    for a in range(n):
        for b in range(a, n):
            M[a, b] = M[b, a] = v[c]; c += 1
    return G(mu, M)
v0 = np.concatenate([mu1, [M1[a, b] for a in range(n) for b in range(a, n)]]); hv = 1e-3; p = len(v0)
Hess = np.zeros((p, p, n))
for i in range(p):
    for j in range(i, p):
        Hess[i, j] = Hess[j, i] = fd(Gvec, v0, i, j, h=hv)
delta_gram = -np.einsum("ijk,ij->k", Hess, Gamma)
print(f"   delta(0) by finite differences in y : {np.array2string(delta_fd, precision=5)}")
print(f"   -<Hess G, Gamma> (theorem)          : {np.array2string(delta_gram, precision=5)}")
print(f"   max abs diff {np.max(np.abs(delta_fd - delta_gram)):.1e}  (|delta| {np.sqrt(np.mean(delta_fd**2)):.2e})")
# 4. Scale-mixture lemma: closure error of the mean vs the kappa4 (2,2) readout, 1-D, Gauss-Hermite in the scale
xg, wg = np.polynomial.hermite_e.hermegauss(60); wg = wg / wg.sum()
for alpha0 in (-1.0, 0.0, 0.7, 1.5):
    sigma0, v = 1.7, 0.01; m0 = alpha0 * sigma0
    def mean_relu(mm, ss): a_ = mm / ss; return ss * (phi(a_) + a_ * ndtr(a_))
    true = np.sum(wg * mean_relu(m0, (1 + np.sqrt(v) * xg) * sigma0))
    s_cl = sigma0 * np.sqrt(1 + v); a_cl = m0 / s_cl; closure = mean_relu(m0, s_cl)
    k4 = 12 * v * sigma0 ** 4; readout = k4 * (a_cl ** 2 - 1) * phi(a_cl) / (24 * s_cl ** 3)
    pred = -0.5 * v * sigma0 * (alpha0 ** 2 - 1) * phi(alpha0)
    print(f"4. scale mixture alpha={alpha0:+.1f}, v={v}: closure-true {closure-true:+.3e}, lemma {pred:+.3e}, with kappa4 readout {closure+readout-true:+.1e} (O(v^2) = {v**2*sigma0:.0e})")
