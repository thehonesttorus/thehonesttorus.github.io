"""Checks for the stage-13 note (v56 read through Klartag-Lehec, Bizeul-Klartag-Lehec, Klartag-Ordentlich, Klartag's
evolving ellipsoid, Sahasrabudhe's survey and the caustic-collar billiard paper).

  1. Eldan's cumulant hierarchy (Bizeul-Klartag-Lehec, Lemma 4.3): along the normalized localization
     dc = A^-1 a dt + A^-1/2 dB, dQ = A^-1 dt, the drift of kappa_m is -(m kappa_m + L_m), with L_3 = 0 and
     L_4(h1..h4) = sum_{j=2..4} <kappa_3(h1,hj,.), A^-1 kappa_3(rest,.)>. Checked on a discrete measure in R^3 by
     Gauss-Hermite averaging over one small time step (two step sizes, Richardson).
  2. Scale-mixture cancellation: for z = (1+d)(mu + s xi) with E d = 0, Var d = eps, the Gaussian closure with the
     mixture's own variance, skewness and kurtosis reproduces E relu z = s m1(mu/s) to O(eps^2).
  3. The dilation sector at width 1024: V_l/4 from the stage-12 increment law, against note XLIV's tau.
  4. Tropical mean width: relu(h_A - h_B) = h_conv(A u B) - h_B, and E f = (E|x|/2)(w(P+) - w(P-)) in R^2.
  5. Zonoid curvature is wall density: Hess_w E relu(w.X) = E[delta(w.X) X X^T] for a non-Gaussian X in R^2.
  6. Negentropy as an integrated localization deficit: D(mu||gamma) = (1/2) int_0^inf (n/(1+t) - E Tr A_t) dt (isotropic mu;
     logistic and two-bump laws in 1D, Monte Carlo over the localization path, t <= 3000, tail O(1/T)).
"""
import itertools
import numpy as np
from math import pi, sqrt, acos, sin, cos, log
from scipy.stats import norm
from scipy.spatial import ConvexHull
from numpy.polynomial.hermite_e import hermegauss

phi, Phi = norm.pdf, norm.cdf
m1 = lambda r: phi(r) + r * Phi(r)

print("== 1. Eldan cumulant hierarchy (drift of kappa_3 and kappa_4)")
rng = np.random.default_rng(3)
d, K = 3, 7
X = rng.normal(size=(K, d)) @ rng.normal(size=(d, d)) + rng.normal(size=(K, d)) ** 2
w0 = rng.random(K) + 0.2
w0 /= w0.sum()


def cumulants(w):
    a = w @ X
    Y = X - a
    A = np.einsum("k,ki,kj->ij", w, Y, Y)
    k3 = np.einsum("k,ki,kj,kl->ijl", w, Y, Y, Y)
    m4 = np.einsum("k,ki,kj,kl,km->ijlm", w, Y, Y, Y, Y)
    k4 = m4 - (np.einsum("ij,lm->ijlm", A, A) + np.einsum("il,jm->ijlm", A, A) + np.einsum("im,jl->ijlm", A, A))
    return a, A, k3, k4


a0, A0, k30, k40 = cumulants(w0)
Ai = np.linalg.inv(A0)
ev, U = np.linalg.eigh(A0)
Aih = U @ np.diag(ev ** -0.5) @ U.T
gx, gw = hermegauss(12)
gw = gw / gw.sum()
nodes = np.array(list(itertools.product(range(12), repeat=d)))
xi = gx[nodes]
wt = np.prod(gw[nodes], axis=1)


def mean_step(dt):
    c = (Ai @ a0) * dt + (xi @ Aih.T) * sqrt(dt)
    Q = Ai * dt
    logw = np.log(w0)[None, :] + c @ X.T - 0.5 * np.einsum("ki,ij,kj->k", X, Q, X)[None, :]
    W = np.exp(logw - logw.max(axis=1, keepdims=True))
    W /= W.sum(axis=1, keepdims=True)
    E3 = np.zeros_like(k30)
    E4 = np.zeros_like(k40)
    for q in range(len(wt)):
        _, _, k3, k4 = cumulants(W[q])
        E3 += wt[q] * k3
        E4 += wt[q] * k4
    return E3, E4


v = np.einsum("ijl,lm->ijm", k30, Ai)
L4 = (np.einsum("ijq,qlm->ijlm", v, k30) + np.einsum("ilq,qjm->ijlm", v, k30) + np.einsum("imq,qjl->ijlm", v, k30))
dts = (2e-4, 1e-4)
D3, D4 = [], []
for dt in dts:
    E3, E4 = mean_step(dt)
    D3.append((E3 - k30) / dt)
    D4.append((E4 - k40) / dt)
R3 = 2 * D3[1] - D3[0]
R4 = 2 * D4[1] - D4[0]
p3, p4 = -3 * k30, -(4 * k40 + L4)
print("  kappa_3 drift: rel. error against -3 kappa_3 = %.2e" % (np.linalg.norm(R3 - p3) / np.linalg.norm(p3)))
print("  kappa_4 drift: rel. error against -(4 kappa_4 + L_4) = %.2e ; |L_4|/|4 kappa_4| = %.3f ; without L_4 the error is %.2e"
      % (np.linalg.norm(R4 - p4) / np.linalg.norm(p4), np.linalg.norm(L4) / np.linalg.norm(4 * k40),
         np.linalg.norm(R4 + 4 * k40) / np.linalg.norm(p4)))
u = rng.normal(size=d)
print("  L_4(u,u,u,u) = %.6f ; 3 <kappa3(u,u,.), A^-1 kappa3(u,u,.)> = %.6f"
      % (np.einsum("ijlm,i,j,l,m->", L4, u, u, u, u), 3 * np.einsum("ijq,qr,lmr,i,j,l,m->", k30, Ai, k30, u, u, u, u)))

print("== 2. Scale-mixture cancellation (two-point and uniform scale laws)")
for law in ("two-point", "uniform"):
    for mu, s in [(0.7, 1.0), (-0.4, 0.8), (1.5, 0.6)]:
        errs = []
        for eps in (0.02, 0.01, 0.005):
            if law == "two-point":
                dv, dw = np.array([-sqrt(eps), sqrt(eps)]), np.array([0.5, 0.5])
            else:
                g, gwt = np.polynomial.legendre.leggauss(20)
                dv, dw = sqrt(3 * eps) * g, gwt / 2
            zx, zw = hermegauss(40)
            zw = zw / zw.sum()
            Z = (1 + dv[:, None]) * (mu + s * zx[None, :])
            Wt = dw[:, None] * zw[None, :]
            Ez = (Wt * Z).sum()
            c = Z - Ez
            k2 = (Wt * c ** 2).sum()
            k3 = (Wt * c ** 3).sum()
            k4 = (Wt * c ** 4).sum() - 3 * k2 ** 2
            S = sqrt(k2)
            R = Ez / S
            g1, g2 = k3 / S ** 3, k4 / S ** 4
            closure = S * (m1(R) - g1 * R * phi(R) / 6 + g2 * (R * R - 1) * phi(R) / 24)
            exact = s * m1(mu / s)  # E(1+d) relu(mu + s xi), by homogeneity (1 + d > 0)
            errs.append(closure - exact)
        print("  %-9s mu=%5.2f s=%.1f: error at eps=.02,.01,.005: %.2e %.2e %.2e (ratios %.2f, %.2f; O(eps^2) -> 4)"
              % (law, mu, s, *errs, errs[0] / errs[1], errs[1] / errs[2]))

print("== 3. The dilation sector at n = 1024")
F = lambda th: acos((sin(th) + (pi - th) * cos(th)) / pi)
J2 = lambda t: 3 * sin(t) * cos(t) + (pi - t) * (1 + 2 * cos(t) ** 2)
th = [pi / 2]
for _ in range(20):
    th.append(F(th[-1]))
n = 1024
nV = 2.0
row = []
for l in range(1, 17):
    nV += 4 * (1.5 - J2(th[l - 1]) / (2 * pi))
    row.append((l, nV / (4 * n)))
print("  V_l/4 (infinite-width increments): " + ", ".join("l=%d %.4f" % r for r in row if r[0] in (4, 6, 8, 10, 12, 14, 16)))
print("  input-radius part 2/(4n) = %.2e ; per-layer increment theta^2/n at l=6,10,14: %.1e %.1e %.1e"
      % (2 / (4 * n), th[5] ** 2 / n, th[9] ** 2 / n, th[13] ** 2 / n))

print("== 4. Tropical mean width in R^2")
rng = np.random.default_rng(11)
k = 8
Wr = rng.normal(size=(k, 2))
vr = rng.normal(size=k)


def zonotope(gens):
    pts = [np.zeros(2)]
    for g in gens:
        pts = pts + [p + g for p in pts]
    return np.array(pts)


Apts = zonotope([vr[i] * Wr[i] for i in range(k) if vr[i] > 0])
Bpts = zonotope([-vr[i] * Wr[i] for i in range(k) if vr[i] < 0])
Cpts = np.vstack([Apts, Bpts])
hsup = lambda P, x: (P @ x.T).max(axis=0)
phis = np.linspace(0, 2 * pi, 20001)[:-1]
U2 = np.stack([np.cos(phis), np.sin(phis)], 1)
f = np.maximum(np.maximum(U2 @ Wr.T, 0) @ vr, 0)
print("  max |relu(h_A - h_B) - (h_conv(AuB) - h_B)| on the circle = %.2e"
      % np.abs(f - (hsup(Cpts, U2) - hsup(Bpts, U2))).max())
per = lambda P: ConvexHull(P).area if len(P) > 2 else 2 * np.linalg.norm(P.max(0) - P.min(0))
Ef_quad = sqrt(pi / 2) * f.mean()
Ef_width = sqrt(pi / 2) / 2 * (per(Cpts) - per(Bpts)) / pi
print("  E f by quadrature %.8f ; (E|x|/2)(w(P+) - w(P-)) with w = perimeter/pi: %.8f" % (Ef_quad, Ef_width))

print("== 5. Zonoid curvature = wall density (2D Gaussian mixture)")
pis = np.array([0.5, 0.3, 0.2])
mus = np.array([[0.8, -0.3], [-1.0, 0.5], [0.2, 1.2]])
Ss = [np.array([[1.0, 0.3], [0.3, 0.6]]), np.array([[0.5, -0.2], [-0.2, 0.9]]), np.array([[0.3, 0.0], [0.0, 0.4]])]


def hfun(w):
    return sum(p * sqrt(w @ S @ w) * m1((w @ m) / sqrt(w @ S @ w)) for p, m, S in zip(pis, mus, Ss))


w = np.array([0.6, -0.9])
e = 1e-4
H = np.zeros((2, 2))
for i in range(2):
    for j in range(2):
        ei, ej = np.eye(2)[i] * e, np.eye(2)[j] * e
        H[i, j] = (hfun(w + ei + ej) - hfun(w + ei - ej) - hfun(w - ei + ej) + hfun(w - ei - ej)) / (4 * e * e)
perp = np.array([-w[1], w[0]]) / np.linalg.norm(w)
ts = np.linspace(-12, 12, 40001)
pts = ts[:, None] * perp[None, :]
dens = sum(p * np.exp(-0.5 * np.einsum("ki,ij,kj->k", pts - m, np.linalg.inv(S), pts - m)) / (2 * pi * sqrt(np.linalg.det(S)))
           for p, m, S in zip(pis, mus, Ss))
Wd = np.einsum("k,ki,kj->ij", dens, pts, pts) * (ts[1] - ts[0]) / np.linalg.norm(w)
print("  Hessian by finite differences:", np.round(H.ravel(), 6), " wall integral:", np.round(Wd.ravel(), 6))

print("== 6. Negentropy as the integrated localization deficit (1D, isotropic)")
rng = np.random.default_rng(2)
for name in ("logistic", "two-bump"):
    if name == "logistic":
        xs = np.linspace(-30, 30, 12001)
        px = 1 / np.cosh(xs / 2) ** 2
    else:
        xs = np.linspace(-6, 6, 6001)
        px = np.exp(-(xs - 1) ** 2 / (2 * 0.25)) + np.exp(-(xs + 1) ** 2 / (2 * 0.25))
    dx = xs[1] - xs[0]
    px /= px.sum() * dx
    mean = (xs * px).sum() * dx
    sd = sqrt(((xs - mean) ** 2 * px).sum() * dx)
    xs = (xs - mean) / sd
    px = px * sd
    dx = xs[1] - xs[0]
    ent = -(px[px > 0] * np.log(px[px > 0])).sum() * dx
    D = 0.5 * log(2 * pi * np.e) - ent
    tgrid = np.concatenate([np.linspace(0, 1, 41)[1:], np.geomspace(1.05, 3000, 160)])
    sel = np.abs(xs) < 12
    xs, px = xs[sel], px[sel]
    dx = xs[1] - xs[0]
    defs = []
    cdf = np.cumsum(px) * dx
    for t in tgrid:
        Xs = np.interp(rng.random(4000), cdf / cdf[-1], xs)
        th_ = t * Xs + sqrt(t) * rng.normal(size=Xs.size)
        lw = th_[:, None] * xs[None, :] - t * xs[None, :] ** 2 / 2 + np.log(np.maximum(px, 1e-300))[None, :]
        Wl = np.exp(lw - lw.max(1, keepdims=True))
        Wl /= Wl.sum(1, keepdims=True)
        mA = Wl @ xs
        vA = Wl @ xs ** 2 - mA ** 2
        defs.append(1 / (1 + t) - vA.mean())
    defs = np.array(defs)
    integ = np.trapezoid(np.concatenate([[0.0], defs]), np.concatenate([[0.0], tgrid]))
    print("  %-8s D(mu||gamma) = %.4f ; (1/2) int (1/(1+t) - E A_t) dt = %.4f (Monte Carlo over theta, t <= 3000)"
          % (name, D, 0.5 * integ))
