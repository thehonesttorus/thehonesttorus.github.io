"""Exactness checks for the identities of DESIGN.md (toy scale).

C1  E2 (barycentre = facet integral) for a 2-D cone, against quadrature.
C2  E f = E[Delta f] for a width-8 depth-2 network: Delta f is a sum of facet terms; checked by smoothing
    ReLU -> softplus_T and Monte Carlo of the Laplacian, T -> 0.
C3  Radial factorisation: a_l(x) = R * a_l(sqrt(n) theta), R = |x|/sqrt(n) independent of theta;
    E a = E R * E_S b, Cov(a) = Cov_S(b) + Var(R) E_S b E_S b^T   (E R^2 = 1).
C4  E3/E4 algebra: Cov(a) = D_m Gamma D_m + D_m X + X^T D_m + Xi on samples (exact identity).
C5  Edgeworth ReLU readout (facet-density form) against numerical integration of an Edgeworth density.
"""
import numpy as np
from scipy import integrate, special

rng = np.random.default_rng(0)
ok = True


def report(name, err, tol):
    global ok
    flag = "OK " if err < tol else "BAD"
    ok &= err < tol
    print(f"[{flag}] {name}: error {err:.2e} (tol {tol:.0e})")


# C1: cone F = {x: <x,u1> > 0, <x,u2> > 0} in R^2; E[x 1_F] = -sum_e nu_e gamma_1(e).
th = 1.1
u1 = np.array([1.0, 0.0]); u2 = np.array([np.cos(th), np.sin(th)])
# the cone between the rays perpendicular to u1, u2; quadrature in polar coordinates
f = lambda t: np.array([np.cos(t), np.sin(t)])
ins = lambda t: (f(t) @ u1 > 0) & (f(t) @ u2 > 0)
ts = (np.arange(4_000_000) + 0.5) * (2 * np.pi / 4_000_000)
dirs = np.stack([np.cos(ts), np.sin(ts)], 1)
tt = ts[(dirs @ u1 > 0) & (dirs @ u2 > 0)]
Er = np.sqrt(np.pi / 2)            # E|x| restricted density: int_0^inf r * r e^{-r^2/2} dr / (2 pi) = sqrt(pi/2)/(2 pi) per angle
bary = np.array([np.cos(tt).sum(), np.sin(tt).sum()]) * (2 * np.pi / 4_000_000) * Er / (2 * np.pi)
# facets: the rays {x: <x,u_i> = 0} on the cone boundary; outward normal of facet i is -u_i; Gaussian measure of a ray = 1/(2 sqrt(2 pi))
g1 = 1 / (2 * np.sqrt(2 * np.pi))
pred = (u1 + u2) * g1
report("C1 barycentre of a 2-D cone = facet integral", np.abs(bary - pred).max(), 1e-4)

# C2: E f = E[Delta f] via softplus smoothing.
n = 8
W1 = rng.standard_normal((n, n)) * np.sqrt(2 / n); W2 = rng.standard_normal((n, n)) * np.sqrt(2 / n)
N = 2_000_000
x = rng.standard_normal((N, n))
z1 = x @ W1; a1 = np.maximum(z1, 0); z2 = a1 @ W2; a2 = np.maximum(z2, 0)
Ef = a2.mean(0)
for T in [0.1, 0.02]:
    sp = lambda u: T * np.logaddexp(0, u / T)
    s1 = lambda u: special.expit(u / T)
    s2 = lambda u: special.expit(u / T) * (1 - special.expit(u / T)) / T
    h1 = sp(z1); y2 = h1 @ W2
    # Laplacian of softplus(y2_j) with y2 = softplus(xW1) W2:
    # grad y2_j = W1 diag(s1(z1)) W2[:, j];  Lap y2_j = sum_i W2_ij s2(z1_i) |W1[:, i]|^2
    c1 = (W1 ** 2).sum(0)
    lap_y2 = (s2(z1) * c1) @ W2
    G = s1(z1)[:, :, None] * W2[None, :, :]                 # (N, i, j) = d y2_j / d h... chain through W1
    gn2 = np.einsum("nij,nkj,ik->nj", G, G, W1.T @ W1, optimize=True)  # |grad y2_j|^2
    lap = s2(y2) * gn2 + s1(y2) * lap_y2
    xg = s1(y2) * ((s1(z1) * z1) @ W2)                     # x . grad softplus(y2_j)
    se = np.sqrt((xg - lap).var(0).max() / N) * 5
    report(f"C2 Stein: E[x.grad f] = E[Lap f] (softplus T={T})", np.abs(xg.mean(0) - lap.mean(0)).max(), se)
    report(f"C2 T->0: E[x.grad f] -> E f for ReLU (T={T})", np.abs(xg.mean(0) - Ef).max(), 3e-2)
xg0 = ((z1 > 0) * z1) @ W2 * (z2 > 0)
report("C2 Euler x.grad f = f for the ReLU net (samplewise)", np.abs(xg0 - a2).max(), 1e-10)

# C3: radial factorisation, width 16 depth 3.
n = 16; L = 3
Ws = [rng.standard_normal((n, n)) * np.sqrt(2 / n) for _ in range(L)]
def net(h):
    for Wl in Ws:
        h = np.maximum(h @ Wl, 0)
    return h
N = 2_000_000
x = rng.standard_normal((N, n)); r = np.linalg.norm(x, axis=1); u = x * (np.sqrt(n) / r)[:, None]
a = net(x); b = net(u)
R = r / np.sqrt(n)
report("C3 a(x) = R b(theta) samplewise", np.abs(a - R[:, None] * b).max(), 1e-9)
ER = np.sqrt(2 / n) * np.exp(special.gammaln((n + 1) / 2) - special.gammaln(n / 2))
x2 = rng.standard_normal((N, n)); b2 = net(x2 * (np.sqrt(n) / np.linalg.norm(x2, axis=1))[:, None])  # independent sphere sample
se = a.std(0).max() / np.sqrt(N) * 4
report("C3 E a = E R * E_S b (independent samples)", np.abs(a.mean(0) - ER * b2.mean(0)).max(), se)
covpred = np.cov(b2.T) + (1 - ER ** 2) * np.outer(b2.mean(0), b2.mean(0))
report("C3 Cov a = Cov_S b + Var R mu mu^T", np.abs(np.cov(a.T) - covpred).max(), 6 * a.var(0).max() / np.sqrt(N))

# C4: barycentric decomposition algebra on samples.
z = x @ Ws[0] - 0.3
g = (z > 0).astype(float); act = np.maximum(z, 0)
p = g.mean(0); Ea = act.mean(0); m = Ea / p
Gam = np.cov(g.T, bias=True)
rr = act - m * g
X = (g - p).T @ (rr - rr.mean(0)) / N
Xi = np.cov(rr.T, bias=True)
lhs = np.cov(act.T, bias=True)
rhs = np.diag(m) @ Gam @ np.diag(m) + np.diag(m) @ X + X.T @ np.diag(m) + Xi
report("C4 Cov(a) = DmGDm + DmX + X^TDm + Xi", np.abs(lhs - rhs).max(), 1e-10)
report("C4 E[r_k | g_k] = 0  (Cov(g_k, r_k) = 0)", np.abs(np.diag(X)).max(), 1e-10)

# C5: Edgeworth ReLU readout.
def readout(mu, var, k3, k4):
    s = np.sqrt(var); t = mu / s; f = np.exp(-t * t / 2) / np.sqrt(2 * np.pi)
    return mu * special.ndtr(t) + s * f - k3 * t * f / (6 * s ** 2) + k4 * (t * t - 1) * f / (24 * s ** 3) \
        + k3 ** 2 * (t ** 4 - 6 * t * t + 3) * f / (72 * s ** 5)
He = lambda k, u: special.eval_hermitenorm(k, u)
for (mu, var, k3, k4) in [(0.3, 1.7, 0.2, 0.1), (-0.5, 2.0, -0.3, 0.4)]:
    s = np.sqrt(var)
    dens = lambda zz: np.exp(-((zz - mu) / s) ** 2 / 2) / (s * np.sqrt(2 * np.pi)) * (
        1 + k3 / (6 * s ** 3) * He(3, (zz - mu) / s) + k4 / (24 * s ** 4) * He(4, (zz - mu) / s)
        + k3 ** 2 / (72 * s ** 6) * He(6, (zz - mu) / s))
    num = integrate.quad(lambda zz: zz * dens(zz), 0, mu + 30 * s, limit=400)[0]
    report(f"C5 Edgeworth ReLU readout mu={mu}", abs(num - readout(mu, var, k3, k4)), 1e-9)

print("ALL OK" if ok else "SOME CHECKS FAILED")
