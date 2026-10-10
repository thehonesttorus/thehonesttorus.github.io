"""Identity checks for the stage-11 note (tilings, Bratteli diagrams, singular foliations).

Each block prints a few numbers that the note quotes. Nothing here is a measurement of an
estimator; these are checks of exact or asymptotic identities.

  1. Bratteli-backprop: forward activations and backward sensitivities are the two harmonic
     states of the layered weighted diagram, their pairing is the output at every layer, and
     per-neuron rescaling (the diagonal flow) preserves every rectangle h_{l,i} d_{l,i}.
  2. Rauzy-Veech induction is a gated linear map: linear on each Rauzy cylinder, with the
     cylinder partition as its activation tiling.
  3. The folding map is a parabolic germ, theta' = theta - theta^2/(3 pi) - theta^3/(18 pi^2) + ...,
     with Fatou coordinate 3 pi/theta + (3/2) log theta + C + O(theta).
  4. Two-tier collapse: the unresolved gate variance at (depth l, localization time t) depends on
     the single Fatou time tau = Phi(theta_0(t)) + l - 1; asymptotically one layer of depth is worth
     sqrt(2)/(3 pi) of sqrt(1+t).
  5. Wall foliation: at a crossing of a deep wall with a wall that feeds it, the formal transverse
     model is the three-line arrangement; linear isotropy drops from dimension 2 to 1 and the
     isotropy Lie algebra becomes aff(1).
  6. Boundary layer: near a single wall, E f(m + sqrt(hbar) Z) - f(m) = sqrt(hbar) kappa psi(d/sqrt(hbar)),
     psi(u) = phi(u) - |u| Phi(-|u|), up to exponentially small terms.
  7. Mean as a border functional in two dimensions: E f = sum over rays of kappa / (2 sqrt(2 pi)).
"""
import numpy as np
import mpmath as mp
from math import erf, sqrt, pi, exp, log, acos, cos, sin

rng = np.random.default_rng(11)
relu = lambda z: np.maximum(z, 0.0)
Phi_n = lambda u: 0.5 * (1.0 + erf(u / sqrt(2.0)))
phi_n = lambda u: exp(-u * u / 2.0) / sqrt(2.0 * pi)


def he_net(n_in, widths, rng):
    Ws, d = [], n_in
    for m in widths:
        Ws.append(rng.normal(0.0, sqrt(2.0 / d), size=(m, d)))
        d = m
    return Ws


def forward(Ws, x):
    hs, zs = [x], []
    for W in Ws:
        z = W @ hs[-1]
        zs.append(z)
        hs.append(relu(z))
    return hs, zs


print("== 1. Bratteli-backprop: <h_l, delta_l> = f at every layer; rescaling invariance")
n, L = 64, 8
Ws = he_net(n, [n] * L, rng)
c = rng.normal(size=n)
x = rng.normal(size=n)
hs, zs = forward(Ws, x)
f = c @ hs[-1]
deltas = [None] * (L + 1)
deltas[L] = c.copy()
for l in range(L, 0, -1):
    g = (zs[l - 1] > 0).astype(float)
    deltas[l - 1] = Ws[l - 1].T @ (g * deltas[l])
pair = [float(hs[l] @ deltas[l]) for l in range(L + 1)]
print("  f = %.12f ; max_l |<h_l,delta_l> - f| = %.2e" % (f, max(abs(p - f) for p in pair)))
# per-neuron rescaling at layer 3: W_3 -> diag(s) W_3, W_4 -> W_4 diag(1/s)
s = np.exp(rng.normal(size=n))
Ws2 = [W.copy() for W in Ws]
Ws2[2] = s[:, None] * Ws2[2]
Ws2[3] = Ws2[3] / s[None, :]
hs2, zs2 = forward(Ws2, x)
d2 = [None] * (L + 1)
d2[L] = c.copy()
for l in range(L, 0, -1):
    g = (zs2[l - 1] > 0).astype(float)
    d2[l - 1] = Ws2[l - 1].T @ (g * d2[l])
rect = hs[3] * deltas[3]
rect2 = hs2[3] * d2[3]
print("  output change %.2e ; max rectangle change at the rescaled layer %.2e ; widths scaled by s: %s"
      % (abs(c @ hs2[-1] - f), np.abs(rect - rect2).max(), np.allclose(hs2[3], s * hs[3])))

print("== 2. Rauzy-Veech induction as a gated linear map (d = 4, permutation (4 3 2 1))")
A = ["A", "B", "C", "D"]


def rv_step(top, bot, lam):
    a0, a1 = top[-1], bot[-1]
    g = lam[a0] > lam[a1]                   # the gate: sign of a linear form
    lam = dict(lam)
    top, bot = list(top), list(bot)
    if g:   # type 0: winner a0, loser a1
        lam[a0] = lam[a0] - lam[a1]
        k = bot.index(a0)
        bot = bot[:k + 1] + [a1] + [b for b in bot[k + 1:] if b != a1]
    else:   # type 1: winner a1, loser a0
        lam[a1] = lam[a1] - lam[a0]
        k = top.index(a1)
        top = top[:k + 1] + [a0] + [t for t in top[k + 1:] if t != a0]
    return top, bot, lam, int(g)


def rv(lamv, k):
    top, bot = ["A", "B", "C", "D"], ["D", "C", "B", "A"]
    lam = dict(zip(A, lamv))
    gates = []
    for _ in range(k):
        top, bot, lam, g = rv_step(top, bot, lam)
        gates.append(g)
    return np.array([lam[a] for a in A]), tuple(gates)


k = 12
lam0 = rng.random(4)
out0, path0 = rv(lam0, k)
# linearity on the cylinder: Jacobian by finite differences, then test at a nearby point of the same cylinder
eps = 1e-7
J = np.zeros((4, 4))
for j in range(4):
    e = np.zeros(4)
    e[j] = eps
    J[:, j] = (rv(lam0 + e, k)[0] - out0) / eps
lin_err, same = [], 0
for _ in range(200):
    lam1 = lam0 + 1e-3 * rng.normal(size=4)
    out1, path1 = rv(lam1, k)
    if path1 == path0:
        same += 1
        lin_err.append(np.abs(out1 - (out0 + J @ (lam1 - lam0))).max())
Jr = np.round(J).astype(int)
print("  Rauzy path of length %d: %s ; Jacobian is an integer matrix: %s ; det = %d"
      % (k, "".join(map(str, path0)), np.allclose(J, Jr, atol=1e-5), round(np.linalg.det(Jr))))
print("  %d/200 perturbed points stay in the cylinder; max deviation from the linear map %.1e"
      % (same, max(lin_err)))
cyl = {}
for _ in range(20000):
    lv = rng.random(4)
    p = rv(lv, 6)[1]
    cyl[p] = cyl.get(p, 0) + 1
print("  distinct gate patterns (cylinders) seen at depth 6: %d of 2^6 = 64" % len(cyl))

print("== 3. Folding map as a parabolic germ; Fatou coordinate")
mp.mp.dps = 40


def Fmap(th):
    one_m = (mp.pi * 2 * mp.sin(th / 2) ** 2 - (mp.sin(th) - th * mp.cos(th))) / mp.pi
    return 2 * mp.asin(mp.sqrt(one_m / 2))


for th in [mp.mpf("0.01"), mp.mpf("0.001")]:
    ser = th - th ** 2 / (3 * mp.pi) - th ** 3 / (18 * mp.pi ** 2)
    print("  theta=%s: (F - series)/theta^4 = %s" % (th, mp.nstr((Fmap(th) - ser) / th ** 4, 6)))
PhiA = lambda th: 3 * mp.pi / th + mp.mpf(3) / 2 * mp.log(th)
th = mp.pi / 2
for l in range(1, 100001):
    th = Fmap(th)
    if l in (10, 100, 1000, 10000, 100000):
        print("  l=%6d theta_l=%s  Phi(theta_l)-l=%s  (3pi/theta_l - l)=%s  3pi/l-law error=%s"
              % (l, mp.nstr(th, 6), mp.nstr(PhiA(th) - l, 8), mp.nstr(3 * mp.pi / th - l, 6),
                 mp.nstr(abs(3 * mp.pi / l - th) / th, 3)))


def Phi_exact(th, N=20000):
    t = mp.mpf(th)
    for _ in range(N):
        t = Fmap(t)
    return PhiA(t) - N


print("== 4. Two-tier collapse in the Fatou time")
mp.mp.dps = 25
# stage-10 table (width 512, 64 pairs): measured / predicted disagreement fractions
layers = [1, 2, 4, 8, 12, 16]
table = {1: [(.333, .333), (.287, .292), (.223, .236), (.174, .173), (.130, .138), (.136, .115)],
         9: [(.140, .144), (.135, .136), (.116, .124), (.106, .105), (.082, .092), (.098, .081)],
         99: [(.044, .045), (.042, .044), (.041, .043), (.040, .041), (.033, .039), (.043, .037)]}
rows = []
for t, vals in table.items():
    th0 = mp.acos(mp.mpf(t) / (1 + t))
    P0 = Phi_exact(th0, 6000)
    for l, (meas, pred) in zip(layers, vals):
        tau = P0 + l - 1
        th = th0
        for _ in range(l - 1):
            th = Fmap(th)
        rows.append((float(tau), t, l, meas, float(th / mp.pi)))
rows.sort()
print("  tau      t   l   measured  predicted theta_{l-1}/pi   (sorted by Fatou time)")
for r in rows:
    print("  %6.2f  %3d  %2d   %.3f     %.3f" % r)
# trade rate: sqrt(1+t_eff) - sqrt(1+t) per layer at large t
for t in [99, 999, 9999]:
    th0 = mp.acos(mp.mpf(t) / (1 + t))
    th = th0
    for _ in range(15):
        th = Fmap(th)
    teff = 1 / (1 - mp.cos(th)) - 1
    rate = (mp.sqrt(1 + teff) - mp.sqrt(1 + t)) / 15
    print("  t=%5d: layer 16 resolved like layer 1 at t_eff=%s; d sqrt(1+t)/d layer = %s (sqrt2/3pi = %s)"
          % (t, mp.nstr(teff, 6), mp.nstr(rate, 5), mp.nstr(mp.sqrt(2) / (3 * mp.pi), 5)))

print("== 5. Wall foliation: transverse isotropy at normal versus bent crossings")


def tangent_space_dim(rays, deg):
    # homogeneous polynomial vector fields of degree deg in R^2, tangent to every ray
    mons = [(i, deg - i) for i in range(deg + 1)]
    rows_ = []
    for r in rays:
        nr = np.array([-r[1], r[0]])
        for tt in np.linspace(0.3, 2.0, 2 * deg + 3):
            p = tt * np.array(r)
            mv = np.array([p[0] ** a * p[1] ** b for a, b in mons])
            rows_.append(np.concatenate([nr[0] * mv, nr[1] * mv]))
    M = np.array(rows_)
    return 2 * len(mons) - np.linalg.matrix_rank(M, tol=1e-9)


for w in [0.0, 0.3]:
    rays = [(0, 1), (0, -1), (-1, 0), (1 / sqrt(1 + w * w), -w / sqrt(1 + w * w))]
    print("  bending w=%.1f: dim of tangent fields of degree 1 = %d, degree 2 = %d, degree 3 = %d"
          % (w, tangent_space_dim(rays, 1), tangent_space_dim(rays, 2), tangent_space_dim(rays, 3)))
# Lie algebra at the bent crossing: E = Euler, theta of degree 2 with det(E, theta) = y1 y2 (y2 + w y1)
import sympy as sp
y1, y2, w_ = sp.symbols("y1 y2 w")
E = sp.Matrix([y1, y2])
# theta = (a(y), b(y)) homogeneous quadratic with y1*b - y2*a = y1*y2*(y2+w*y1); choose a = -y2*(y2+w*y1)... solve generically
a0, a1, a2, b0, b1, b2 = sp.symbols("a0 a1 a2 b0 b1 b2")
ath = a0 * y1 ** 2 + a1 * y1 * y2 + a2 * y2 ** 2
bth = b0 * y1 ** 2 + b1 * y1 * y2 + b2 * y2 ** 2
eqs = sp.Poly(sp.expand(y1 * bth - y2 * ath - y1 * y2 * (y2 + w_ * y1)), y1, y2).coeffs()
sol = sp.solve(eqs, [a0, a1, a2, b0, b1, b2], dict=True)[0]
th_ = sp.Matrix([ath, bth]).subs(sol).subs({s_: 0 for s_ in [a0, a1, a2, b0, b1, b2]})


def bracket(X, Y):
    JX = X.jacobian([y1, y2])
    JY = Y.jacobian([y1, y2])
    return sp.simplify(JY * X - JX * Y)


br = bracket(E, th_)
print("  theta =", list(th_), "; [E, theta] - theta =", list(sp.simplify(br - th_)),
      "; tangent to the phantom line y2=0, y1>0:", sp.simplify(th_[1].subs(y2, 0)) == 0)

print("== 6. Boundary layer of the heat deformation near one wall")
n2 = 2
ang = np.arange(6) * pi / 6 + rng.uniform(-0.1, 0.1, size=6)
W1 = np.stack([np.cos(ang), np.sin(ang)], 1) * rng.uniform(0.5, 1.5, size=(6, 1))
w2 = rng.normal(0, 1, size=6)
fnet = lambda X: relu(X @ W1.T) @ w2       # X: (N,2)
# choose a point near the wall of unit 0, away from the other walls
u0 = W1[0] / np.linalg.norm(W1[0])
tang = np.array([-u0[1], u0[0]])
best = None
for r in np.linspace(1.0, 8.0, 141):
    for sgn in (1, -1):
        p = sgn * r * tang
        dist_other = min(abs(W1[j] @ p) / np.linalg.norm(W1[j]) for j in range(1, 6))
        if best is None or dist_other > best[0]:
            best = (dist_other, p)
dist_other, p0 = best
kappa = w2[0] * np.linalg.norm(W1[0])
psi = lambda u: phi_n(u) - abs(u) * Phi_n(-abs(u))
gh_x, gh_w = np.polynomial.hermite_e.hermegauss(200)
gh_w = gh_w / gh_w.sum()
print("  distance from the base point to the other walls: %.2f" % dist_other)
for hbar in [0.01, 0.05, 0.2, 0.5]:
    sd = sqrt(hbar)
    for dsig in [0.0, 0.7]:
        m = p0 + dsig * sd * u0
        # 2D quadrature along the normal (exact kink handled by splitting) and tangential Gauss-Hermite
        tot = 0.0
        for tx, tw in zip(gh_x, gh_w):
            # normal direction: integrate f(m + sd*(s*u0 + tx*tang)) against N(0,1) in s with a fine grid
            s = np.linspace(-9, 9, 36001)
            pts = m[None, :] + sd * (s[:, None] * u0[None, :] + tx * tang[None, :])
            vals = fnet(pts) * np.exp(-s ** 2 / 2) / sqrt(2 * pi)
            tot += tw * np.trapezoid(vals, s)
        lhs = tot - fnet(m[None, :])[0]
        rhs = sd * kappa * psi(dsig)
        print("  hbar=%.2f d/sqrt(hbar)=%.1f: F-f=%.6e  sqrt(hbar) kappa psi=%.6e  rel.err=%.1e"
              % (hbar, dsig, lhs, rhs, abs(lhs - rhs) / abs(rhs)))

print("== 7. Mean as a border functional (two-dimensional input, rays as walls)")
W1b = rng.normal(0, 1, size=(16, 2))
W2b = rng.normal(0, sqrt(2 / 16), size=(16, 16))
cb = rng.normal(size=16)
fb = lambda X: relu(relu(X @ W1b.T) @ W2b.T) @ cb
th = np.linspace(0, 2 * pi, 2_000_001)
U = np.stack([np.cos(th), np.sin(th)], 1)
fv = fb(U)
mean_exact = sqrt(pi / 2) * np.trapezoid(fv, th) / (2 * pi)
# exact breakpoints: sign changes of every pre-activation along the circle, refined by bisection;
# kappa = jump of the normal derivative, from the backprop gradient on either side


def grad_b(X):
    z1 = X @ W1b.T
    h1 = relu(z1)
    z2 = h1 @ W2b.T
    d2 = cb * (z2 > 0)
    d1 = (d2 @ W2b) * (z1 > 0)
    return d1 @ W1b


def preacts(t):
    X = np.array([[cos(t), sin(t)]])
    z1 = X @ W1b.T
    return np.concatenate([z1[0], (relu(z1) @ W2b.T)[0]])


Z = np.concatenate([U @ W1b.T, relu(U @ W1b.T) @ W2b.T], 1)
kap = []
for j in range(Z.shape[1]):
    for i in np.where(np.sign(Z[:-1, j]) * np.sign(Z[1:, j]) < 0)[0]:
        a_, b_ = th[i], th[i + 1]
        for _ in range(60):
            mid = 0.5 * (a_ + b_)
            if np.sign(preacts(mid)[j]) == np.sign(preacts(a_)[j]):
                a_ = mid
            else:
                b_ = mid
        ts = 0.5 * (a_ + b_)
        e_t = np.array([-sin(ts), cos(ts)])
        gp = grad_b(np.array([[cos(ts + 1e-9), sin(ts + 1e-9)]]))[0]
        gm = grad_b(np.array([[cos(ts - 1e-9), sin(ts - 1e-9)]]))[0]
        kap.append((gp - gm) @ e_t)
mean_border = sum(kap) / (2 * sqrt(2 * pi))
print("  breakpoints on the circle: %d ; E f by angular quadrature = %.8f ; sum_rays kappa/(2 sqrt(2pi)) = %.8f"
      % (len(kap), mean_exact, mean_border))
