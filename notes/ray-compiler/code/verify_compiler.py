# Checks of the supplied context-typed compiler and suffix-risk results (note: notes/ray-compiler/README.md).
#   python verify_compiler.py            (seconds; exact rational arithmetic and low-dimensional quadrature only)
# V1 terminal cosine-transform multipliers of relu on S^(n-1); V2 the radial moment recurrence and the two-scalar reduction;
# V3 the degree-p null tuples, as ray integrals; V4 the inner-context coefficient; V5 the minimum-score gauge equation;
# V6 the finite-width two-input transition at c = 1 and its Cauchy-Schwarz contraction.
import itertools, math
from fractions import Fraction as Fr
import numpy as np
from scipy import integrate

rng = np.random.default_rng(7)
out = []
say = lambda s: (print(s), out.append(s))

# ---------- V1: Funk-Hecke multipliers mu_l = E[relu(t) P_l(t)], t the first coordinate of a uniform point of S^(n-1) ----------
def multipliers(n, L):
    # normalised Gegenbauer P_l (P_l(1) = 1): (l + n - 2) P_(l+1) = (2l + n - 2) t P_l - l P_(l-1); coefficients as Fractions
    P = [[Fr(1)], [Fr(0), Fr(1)]]
    for l in range(1, L):
        a, b = Fr(2 * l + n - 2, l + n - 2), Fr(l, l + n - 2)
        nxt = [Fr(0)] * (l + 2)
        for j, c in enumerate(P[l]):
            nxt[j + 1] += a * c
        for j, c in enumerate(P[l - 1]):
            nxt[j] -= b * c
        P.append(nxt)
    ev = lambda m: math.prod((Fr(2 * i + 1, n + 2 * i) for i in range(m)), start=Fr(1))       # E t^(2m)
    ab = lambda m: math.prod((Fr(2 * i + 2, n + 2 * i + 1) for i in range(m)), start=Fr(1))   # E|t| t^(2m) / E|t|
    mu = []                    # (rational part, coefficient of E|t|):  E relu(t) t^j = E t^(j+1)/2 + E|t| t^j / 2
    for l in range(L + 1):
        lin = sum((c * ev((j + 1) // 2) / 2 for j, c in enumerate(P[l]) if j % 2 == 1), Fr(0))
        absp = sum((c * ab(j // 2) / 2 for j, c in enumerate(P[l]) if j % 2 == 0), Fr(0))
        mu.append((lin, absp))
    return mu

def ncount(n, l):
    return math.comb(n + l - 1, l) - (math.comb(n + l - 3, l - 2) if l >= 2 else 0)

say("V1 terminal multipliers of relu on S^(n-1): mu_1 = 1/(2n), odd l >= 3 vanish, mu_(2k+2)/mu_(2k) = -(2k-1)/(n+2k+1)")
for n in (3, 8, 64, 1024):
    mu = multipliers(n, 12)
    ok1 = mu[1] == (Fr(1, 2 * n), Fr(0))
    odd = all(mu[l] == (Fr(0), Fr(0)) for l in range(3, 13, 2))
    rat = [mu[2 * k + 2][1] / mu[2 * k][1] for k in range(1, 6)]
    law = all(rat[k - 1] == Fr(-(2 * k - 1), n + 2 * k + 1) for k in range(1, 6))
    e4 = float((mu[4][1] / mu[2][1]) ** 2 * ncount(n, 4) / ncount(n, 2))
    say(f"  n={n:5d}: mu_1 = 1/(2n) {ok1}; odd l=3..11 zero {odd}; ratio law k=1..5 {law}; mu_4/mu_2 = {rat[0]}; "
        f"energy ratio (mu_4^2 N_4)/(mu_2^2 N_2) = {e4:.4f}")
say("  (the energy a single unit puts in degree 4 is ~1/12 of degree 2 at every n: the 1/(n+3) of the multiplier is the\n"
    "   harmonic normalisation, not a suppression of fourth-order directional structure)")

# ---------- V2: radial moments m_k = int_0^inf r^(n+k) exp(-r^2/2 + beta r) dr ----------
def mk(n, k, beta):
    f = lambda r: math.exp((n + k) * math.log(r) - r * r / 2 + beta * r - ((n + k) * math.log(rs) - rs * rs / 2 + beta * rs)) if r > 0 else 0.0
    rs = (beta + math.sqrt(beta * beta + 4 * (n + k))) / 2        # mode, for scaling
    v, _ = integrate.quad(f, 0, rs + 40, limit=400, epsabs=0, epsrel=1e-13)
    return v * math.exp((n + k) * math.log(rs) - rs * rs / 2 + beta * rs)

say("V2 radial recurrence m_(k+2) = beta m_(k+1) + (n+k+1) m_k, and the two-scalar reduction of a ray polynomial")
for n, beta in ((5, 0.7), (40, -1.3), (200, 2.1)):
    m = [mk(n, k, beta) for k in range(8)]
    err = max(abs(m[k + 2] - beta * m[k + 1] - (n + k + 1) * m[k]) / m[k + 2] for k in range(6))
    a = rng.standard_normal(8)                        # s(r) = sum a_k r^k  ->  q = c0 m_0 + c1 m_1
    c = list(a)
    for k in range(7, 1, -1):                         # r^k = beta r^(k-1) + (n+k-1) r^(k-2) modulo radial divergences
        c[k - 1] += beta * c[k]; c[k - 2] += (n + k - 1) * c[k]; c[k] = 0.0
    q_direct = sum(a[k] * m[k] for k in range(8)); q_two = c[0] * m[0] + c[1] * m[1]
    say(f"  n={n:3d} beta={beta:+.1f}: recurrence rel. error {err:.1e}; two-scalar reduction rel. error "
        f"{abs(q_direct - q_two) / abs(q_direct):.1e}; kappa_n = m_1/m_0 = {m[1] / m[0]:.4f} (sqrt(n+1/2) + beta/2 = {math.sqrt(n + 0.5) + beta / 2:.4f})")

# ---------- V3: null tuples as ray integrals, n = 3, non-centred correlated reference ----------
nd = 3
A = rng.standard_normal((nd, nd)); Sig = A @ A.T / nd + 0.5 * np.eye(nd); Pm = np.linalg.inv(Sig)
mu = np.array([0.6, -0.3, 0.45]); G = rng.standard_normal((nd, nd)); G = (G + G.T) / 2

def sym(T):
    k = T.ndim
    return sum(np.transpose(T, p) for p in itertools.permutations(range(k))) / math.factorial(k)

def sp(a, B):                                         # normalised symmetric product of a tensor and a tensor
    return sym(np.multiply.outer(a, B))

def dens_derivs(z):                                   # phi and its first four derivative tensors at z
    v = Pm @ (z - mu); ph = math.exp(-0.5 * (z - mu) @ v) / math.sqrt((2 * math.pi) ** nd * np.linalg.det(Sig))
    d2 = np.multiply.outer(v, v) - Pm
    d3 = -(np.einsum("a,b,c->abc", v, v, v) - 3 * sym(np.multiply.outer(Pm, v)))
    d4 = (np.einsum("a,b,c,d->abcd", v, v, v, v) - 6 * sym(np.einsum("ab,c,d->abcd", Pm, v, v))
          + 3 * sym(np.multiply.outer(Pm, Pm)))
    return ph, d2 * ph, d3 * ph, d4 * ph

def ray(T2, T3, T4, p, u):                            # q(u) = int_0^inf r^(n-1+p) dp(r u) dr
    def f(r):
        ph, d2, d3, d4 = dens_derivs(r * u)
        dp = np.sum(T2 * d2) / 2 - np.sum(T3 * d3) / 6 + np.sum(T4 * d4) / 24
        return r ** (nd - 1 + p) * dp
    return integrate.quad(f, 0, 40, limit=400, epsabs=1e-15, epsrel=1e-12)[0]

def scale(T2, T3, T4, p, u):                          # size of the individual terms, for the relative error
    return max(abs(ray(T2, 0 * T3, 0 * T4, p, u)), abs(ray(0 * T2, T3, 0 * T4, p, u)), abs(ray(0 * T2, 0 * T3, T4, p, u)))

say("V3 degree-p null tuples (dSigma, dk3, dk4) = (-(p-2)/12 G, mu.G/4, Sigma.G): ray integrals q(u) = 0 for every direction")
S4, S3 = sp(Sig, G), sp(mu, G)
for p in (1, 2, 3):
    rel = []
    for _ in range(4):
        u = rng.standard_normal(nd); u /= np.linalg.norm(u)
        rel.append(abs(ray(-(p - 2) / 12 * G, S3 / 4, S4, p, u)) / scale(-(p - 2) / 12 * G, S3 / 4, S4, p, u))
    ctrl = abs(ray(-(p - 1) / 12 * G, S3 / 4, S4, p, u)) / scale(-(p - 1) / 12 * G, S3 / 4, S4, p, u)
    say(f"  p={p}: max |q|/term over 4 directions {max(rel):.1e}; control with the degree-(p+1) covariance coefficient {ctrl:.1e}")
t = 0.3                                               # gain mode: (t (Sigma + mu mu^T), 6 t mu.Sigma, 12 t Sigma.Sigma)
u = rng.standard_normal(nd); u /= np.linalg.norm(u)
g2, g3, g4 = t * (Sig + np.outer(mu, mu)), 6 * t * sp(mu, Sig), 12 * t * sp(Sig, Sig)
say(f"  gain mode, p=1: |q|/term {abs(ray(g2, g3, g4, 1, u)) / scale(g2, g3, g4, 1, u):.1e}; p=2 (not null): "
    f"{abs(ray(g2, g3, g4, 2, u)) / scale(g2, g3, g4, 2, u):.1e}")

# ---------- V4: inner-context coefficient: a K4 core inside a degree-m context reads as (p-m-2)/12 ----------
# In one dimension at N(m0, s2), take F = relu and the context of a kappa3 source: g = E-derivative of order 3, degree p - 3.
# Check L4[s2 g4] f''' = ((p-3-2)/12) L2[g4] f''' - (1/4) L3[m0 g4] f''' for f = relu (p = 1) by Gauss-Hermite.
m0, s2, g4 = 0.37, 1.3, 0.8
xg, wg = np.polynomial.hermite_e.hermegauss(80); wg = wg / wg.sum()
phi = lambda x: math.exp(-x * x / 2) / math.sqrt(2 * math.pi)
s = math.sqrt(s2)
def Ed(k):                                            # E relu^(k)(Z), Z ~ N(m0, s2), k >= 2: derivatives of the Gaussian density
    a = -m0 / s                                        # relu'' = delta, so E relu^(k) = (-1)^k d^(k-2) phi_Z(0) / dz^(k-2)
    he = np.polynomial.hermite_e.HermiteE.basis(k - 2)(a)
    return (1 / s) ** (k - 1) * he * phi(a) * (1 if k % 2 == 0 else 1) * (-1) ** 0
lhs = s2 * g4 * Ed(7) / 24 * 6                        # L4[Sigma.G] applied to g = f'''/6 (kappa3 context), up to the common 1/6
rhs = ((1 - 3 - 2) / 12) * g4 * Ed(5) / 2 * 6 - 0.25 * m0 * g4 * Ed(6) / 6 * 6
say(f"V4 inner context (kappa3 x kappa4 term, p=1): coefficient (p-3-2)/12 = -1/3: lhs {lhs:+.6e} rhs {rhs:+.6e} "
    f"rel {abs(lhs - rhs) / abs(lhs):.1e}; with the outer -1/12 instead: {abs(lhs - (((1 - 2) / 12) * g4 * Ed(5) / 2 * 6 - 0.25 * m0 * g4 * Ed(6))) / abs(lhs):.1e}")

# ---------- V5: minimum-score gauge: argmin_G |S2 - G/12|^2/2 + |S3 - a.G/4|^2/6 + |S4 - I.G|^2/24 (whitened) ----------
n5 = 5; a5 = rng.standard_normal(n5) * 0.7
S2 = rng.standard_normal((n5, n5)); S2 = (S2 + S2.T) / 2
S3 = sym(rng.standard_normal((n5,) * 3)); S4 = sym(rng.standard_normal((n5,) * 4)); I5 = np.eye(n5)
idx = [(i, j) for i in range(n5) for j in range(i, n5)]
def Gof(x):
    Gm = np.zeros((n5, n5))
    for (i, j), v in zip(idx, x):
        Gm[i, j] = Gm[j, i] = v
    return Gm
def resid(x):
    Gm = Gof(x)
    return np.concatenate([(S2 - Gm / 12).ravel() / math.sqrt(2), (S3 - sp(a5, Gm) / 4).ravel() / math.sqrt(6),
                           (S4 - sp(I5, Gm)).ravel() / math.sqrt(24)])
J = np.stack([resid(np.eye(len(idx))[k]) - resid(np.zeros(len(idx))) for k in range(len(idx))], 1)
xls = np.linalg.lstsq(J, -resid(np.zeros(len(idx))), rcond=None)[0]; Gls = Gof(xls)
D = S2 + np.einsum("a,abc->bc", a5, S3) + np.einsum("aacd->cd", S4)
c5 = 4 * n5 + 18 + 2 * a5 @ a5
lhs5 = c5 * Gls + 4 * np.trace(Gls) * I5 + 2 * (Gls @ np.outer(a5, a5) + np.outer(a5, a5) @ Gls)
say(f"V5 gauge equation cG + 4 tr(G) I + 2(G a a^T + a a^T G) = 24 D, c = 4n + 18 + 2|a|^2, D = S2 + a.S3 + tr12 S4: "
    f"max |lhs - 24 D| / |24 D| = {np.abs(lhs5 - 24 * D).max() / np.abs(24 * D).max():.1e}")
# closed-form O(n^2) solve: on the complement of a, (c + 2|a|^2 terms) act by scalars; solve in the basis {a a^T, a b^T + b a^T, rest}
def gauge_solve(D, a, n):
    c = 4 * n + 18 + 2 * a @ a; aa = a @ a
    # G a a^T + a a^T G is diagonal in the split (a-parallel, a-perp): eigen-coefficients c + 4aa (par-par), c + 2aa (par-perp), c (perp-perp)
    e = a / math.sqrt(aa); Pp = np.outer(e, e); Q = np.eye(n) - Pp
    Dpp, Dpq, Dqq = Pp @ D @ Pp, Pp @ D @ Q + Q @ D @ Pp, Q @ D @ Q
    # trace couples the pp and qq blocks: solve for t = tr(G) self-consistently
    # G_pp = (24 Dpp - 4 t Pp)/(c + 4 aa), G_qq = (24 Dqq - 4 t Q)/c, G_pq = 24 Dpq/(c + 2 aa)
    t = (24 * np.trace(Dpp) / (c + 4 * aa) + 24 * np.trace(Dqq) / c) / (1 + 4 / (c + 4 * aa) + 4 * (n - 1) / c)
    return (24 * Dpp - 4 * t * Pp) / (c + 4 * aa) + (24 * Dqq - 4 * t * Q) / c + 24 * Dpq / (c + 2 * aa)
say(f"  closed-form O(n^2) solve against least squares: max |G - G_ls| = {np.abs(gauge_solve(D, a5, n5) - Gls).max():.1e}")

# ---------- V6: finite-width transition Q_n f(c) = E[a f(c')], a = (2/n) sqrt(S_A S_B), c' = T / sqrt(S_A S_B) ----------
say("V6 finite-width transition: E a = 1 exactly at c = 1 (a = 2 S_A / n); E a < 1 for c < 1 (Cauchy-Schwarz)")
for n6 in (16, 64, 256):
    for c in (1.0, 0.9, 0.3):
        g1 = rng.standard_normal((20000, n6)); g2 = c * g1 + math.sqrt(max(0.0, 1 - c * c)) * rng.standard_normal((20000, n6))
        ra, rb = np.maximum(g1, 0), np.maximum(g2, 0)
        SA, SB, T = (ra * ra).sum(1), (rb * rb).sum(1), (ra * rb).sum(1)
        aw = 2 / n6 * np.sqrt(SA * SB)
        rho = (math.sqrt(1 - c * c) + (math.pi - math.acos(c)) * c) / math.pi
        say(f"  n={n6:4d} c={c:.1f}: E a = {aw.mean():.5f} +- {aw.std() / math.sqrt(len(aw)):.5f}; E[a c'] = {(aw * T / np.sqrt(SA * SB)).mean():.5f} "
            f"vs rho(c) = {rho:.5f}")

open(__file__.replace("code/verify_compiler.py", "outputs/verify_compiler.txt"), "w").write("\n".join(out) + "\n")
