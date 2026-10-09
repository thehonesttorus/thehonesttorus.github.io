"""Identity checks for notes/stage10 (no experiments: each check compares two sides of a proved identity).
  python scripts/verify_stage10.py"""
import numpy as np
from math import erf, sqrt, pi, log, cos, acos
rng = np.random.default_rng(10)
Phi = lambda t: 0.5 * (1 + np.vectorize(erf)(np.asarray(t) / sqrt(2)))
phi = lambda t: np.exp(-np.asarray(t) ** 2 / 2) / sqrt(2 * pi)

def net(Ws, x):
    zs, h = [], x
    for W in Ws:
        z = h @ W.T; zs.append(z); h = np.maximum(z, 0)
    return zs, h

print("1. kink-current formula (Thm 3.1): E h_L vs sum over walls of E[delta(z)|grad z|^2 (dh_L/dh_l) e_a]")
for n, L, N, eps in ((6, 3, 3_000_000, 0.05), (16, 3, 3_000_000, 0.05)):
    Ws = [rng.standard_normal((n, n)) * sqrt(2 / n) for _ in range(L)]
    x = rng.standard_normal((N, n)); zs, hL = net(Ws, x); g = [(z > 0) for z in zs]
    dead = np.mean(np.any(np.stack([~gg.any(1) for gg in g[:-1]]), 0))
    lhs = hL.mean(0); rhs = np.zeros(n)
    for l in range(L):
        for a in range(n):
            s_ = np.abs(zs[l][:, a]) < eps
            r = np.tile(Ws[l][a], (s_.sum(), 1))                       # row a of dz_l/dx, back-propagated
            for k in range(l - 1, -1, -1): r = (r * g[k][s_]) @ Ws[k]
            v = np.zeros((s_.sum(), n)); v[:, a] = 1.0                  # unit kick at h_{l,a}
            for k in range(l + 1, L):
                v = v @ Ws[k].T
                if k < L - 1: v = v * g[k][s_]
            if l < L - 1: v = v * g[L - 1][s_]
            rhs += (np.sum(r ** 2, 1)[:, None] * v).sum(0) / (N * 2 * eps)
    print("   n=%2d L=%d  P(some layer entirely inactive)=%.1e  max |lhs-rhs|/max|lhs| = %.4f" % (n, L, dead, np.max(np.abs(lhs - rhs)) / np.max(np.abs(lhs))))
    print("        lhs", np.round(lhs[:6], 4), "\n        rhs", np.round(rhs[:6], 4))

print("2. Crofton (Thm 1.2): crossings of random walls by a curve = n * angular length / pi")
m, T = 50, 4000
c = np.cumsum(rng.standard_normal((T, m)) * 0.05, axis=0) + rng.standard_normal(m)
ch = c / np.linalg.norm(c, axis=1, keepdims=True); Lam = np.sum(np.linalg.norm(np.diff(ch, axis=0), axis=1))
w = rng.standard_normal((20000, m)); s = np.sign(c @ w.T); cross = np.sum(s[1:] != s[:-1]) / 20000
print("   mean crossings per wall %.4f   angular length / pi %.4f" % (cross, Lam / pi))

print("3. folding (Thm 1.4): angle recursion vs 3*pi/l and the heat-trace flow (Prop 5.6)")
th = pi / 2
for l in range(1, 2001):
    th = acos(min(1.0, (np.sin(th) + (pi - th) * np.cos(th)) / pi))
    if l in (10, 100, 1000, 2000): print("   l=%5d theta=%.5f  3pi/l=%.5f  ratio %.4f  beta*l^2/(9pi^2/2)=%.4f" % (l, th, 3 * pi / l, th * l / (3 * pi), -log(cos(th)) * l * l / (9 * pi * pi / 2)))

print("4. unit angular gain (Prop 1.3): E|dc'|^2/|dc|^2 and E|c'|^2/|c|^2 at sigma_w^2 = 2")
n2, R = 400, 400; cc = rng.standard_normal(n2); v = rng.standard_normal(n2); v -= v @ cc / (cc @ cc) * cc
g1 = g2 = 0.0
for _ in range(R):
    W = rng.standard_normal((n2, n2)) * sqrt(2 / n2); z = W @ cc; g1 += np.sum((z > 0) * (W @ v) ** 2) / (v @ v); g2 += np.sum(np.maximum(z, 0) ** 2) / (cc @ cc)
print("   tangent gain %.4f  norm gain %.4f" % (g1 / R, g2 / R))

print("5. Mehler-KMS (Prop 5.5) and ReLU Hermite coefficients: E relu(X)relu(Y) vs sum c_k^2 rho^k / k!")
from math import factorial
def dfact(k): return 1 if k <= 0 else k * dfact(k - 2)
ck = [1 / sqrt(2 * pi), 0.5] + [((-1) ** ((k - 2) // 2) * dfact(k - 3) / sqrt(2 * pi) if k % 2 == 0 else 0.0) for k in range(2, 200)]
for rho in (0.3, 0.7, 0.95):
    th = acos(rho); exact = (np.sin(th) + (pi - th) * rho) / (2 * pi)
    from math import lgamma, exp
    series = sum(exp(2 * log(abs(ck[k])) + k * log(rho) - lgamma(k + 1)) for k in range(len(ck)) if ck[k] != 0)
    print("   rho=%.2f closed form %.8f  Hermite series %.8f" % (rho, exact, series))

print("6. contact order law (Thm 5.2): a = int_0^1 r^(d+1) dr = 1/(d+2), checked by the tau-integral of stage 9")
for d in (0, 2, 4):
    He = {0: lambda s: 1 + 0 * s, 2: lambda s: s ** 2 - 1, 4: lambda s: s ** 4 - 6 * s ** 2 + 3}[d]
    taus = np.linspace(0, 400, 400001); wts = 0.5 * (1 + taus) ** -1.5
    gh, gw = np.polynomial.hermite_e.hermegauss(80); gw = gw / gw.sum()
    Ed = np.array([np.sum(gw * phi(np.sqrt(t) * gh) * He(np.sqrt(t) * gh)) for t in taus[::400]])
    Ed = np.interp(taus, taus[::400], Ed); a = np.trapezoid(wts * Ed, taus) / (phi(0.0) * He(0.0))
    print("   d=%d  a=%.4f  1/(d+2)=%.4f" % (d, a, 1 / (d + 2)))

print("7. births are contacts (Thm 4.3): kappa_3 of relu outputs / eps^2 vs the hub (kink) formula as eps -> 0 (Sobol QMC)")
from scipy.stats import qmc, norm as _norm
Sb = np.array([[0, 1, 0.5], [1, 0, -0.7], [0.5, -0.7, 0]]); mvec = np.array([0.3, -0.2, 0.5])
xi = _norm.ppf(qmc.Sobol(3, scramble=True, seed=7).random_base2(23).clip(1e-12, 1 - 1e-12))
for eps2 in (0.16, 0.08, 0.04, 0.02):
    S = np.eye(3) + eps2 * Sb; Lc = np.linalg.cholesky(S); zz = mvec + xi @ Lc.T; hh = np.maximum(zz, 0); hc = hh - hh.mean(0)
    k3 = np.mean(hc[:, 0] * hc[:, 1] * hc[:, 2]); sd = np.sqrt(np.diag(S)); al = mvec / sd; Ph, w2 = Phi(al), phi(al) / sd
    pred = sum(Ph[a] * S[a, c] * w2[c] * S[c, b] * Ph[b] for (a, c, b) in ((1, 0, 2), (0, 1, 2), (0, 2, 1)))
    print("   eps=%.2f  kappa_3/eps^2 = %+.5f   hub formula/eps^2 = %+.5f   ratio %.3f" % (eps2, k3 / eps2 ** 2, pred / eps2 ** 2, k3 / pred))

print("8. complementarity (Thm 4.7): ||E_B E_A - J|| = ||C|| and ||C|| ~ 2 sqrt(2/n) for a He layer covariance")
for nn in (64, 256):
    W = rng.standard_normal((nn, nn)) * sqrt(2 / nn); Sx = W @ W.T; ev, U = np.linalg.eigh(Sx)
    Cm = U ** 2 - 1 / nn; normC = np.linalg.norm(Cm, 2)
    # superoperator norm of E_B E_A - J on L^2(tau): restrict to diagonal inputs Y = diag(y), tau(Y)=0 (E_A kills the rest)
    M = (U ** 2).T            # y -> coefficients <u_i, diag(y) u_i>
    Pm = np.eye(nn) - 1 / nn  # remove the trace part
    op = np.linalg.norm(M @ Pm, 2)
    print("   n=%d  ||C||=%.4f  ||E_B E_A - J|| on diagonal=%.4f  2*sqrt(2/n)=%.4f" % (nn, normC, op, 2 * sqrt(2 / nn)))

print("9. frame bound (Thm 4.4): ||U(I-M)Z^T||^2 >= lmin(U^TU) lmin(Z^TZ) (N-k) for random rank-k M")
Nh, kk = 40, 10; Uu = rng.standard_normal((60, Nh)); Zz = rng.standard_normal((70, Nh)); worst = np.inf
lb = np.linalg.eigvalsh(Uu.T @ Uu)[0] * np.linalg.eigvalsh(Zz.T @ Zz)[0] * (Nh - kk)
for _ in range(2000):
    M = rng.standard_normal((Nh, kk)) @ rng.standard_normal((kk, Nh)) / Nh; worst = min(worst, np.linalg.norm(Uu @ (np.eye(Nh) - M) @ Zz.T) ** 2)
print("   min over 2000 random M: %.2f  >=  bound %.2f" % (worst, lb))
