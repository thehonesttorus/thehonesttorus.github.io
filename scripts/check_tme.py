"""Transcription checks for whest/tme.py against per-orthant Gauss-Legendre quadrature of the exact Gaussian
expectations (ReLU is polynomial on each orthant, so the quadrature converges exponentially).

1. Tokens: G, G2, kappa3(relu), kappa4(relu), F1 = E d(u-mu)^2, F2 = E d^2(u-mu)^2 (Stein form) vs 1-D quadrature.
2. Pairs at small correlation c: two-token cov, coincident kappa3 F1 a1 c + F2 a2 c^2/2, K^u_ab formula; residuals ~c^3.
3. Distinct triple: the cherry a1 a1 a2 C_ac C_bc (+ perms) vs 3-D quadrature; residual ~C^3.
4. Register readouts reproduce the represented tensor's (a,a,a) and (a,a,b) entries (sym = sum of the 3 placements).
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from whest.tme import tokens, readouts

relu = lambda z: np.maximum(z, 0)


def nodes1(m, s, q=120, R=9.0):
    """Gauss-Legendre nodes/weights on [m-Rs, 0] and [0, m+Rs] (split at the kink), with the N(m,s^2) density."""
    xg, wg = np.polynomial.legendre.leggauss(q); zs = []; ws = []
    for lo, hi in ((m - R * s, 0.0), (0.0, m + R * s)):
        if hi <= lo: continue
        z = 0.5 * (hi - lo) * xg + 0.5 * (hi + lo); zs.append(z); ws.append(0.5 * (hi - lo) * wg)
    return np.concatenate(zs), np.concatenate(ws)


def E1(f, m, s):
    z, w = nodes1(m, s); dens = np.exp(-0.5 * ((z - m) / s) ** 2) / (s * np.sqrt(2 * np.pi))
    return np.sum(w * dens * f(z))


def Ek(f, m, C, q=70):
    """E f(z), z ~ N(m, C) in k = 2 or 3 dimensions, box [m_i - 8 s_i, m_i + 8 s_i] split at every z_i = 0."""
    k = len(m); s = np.sqrt(np.diag(C)); grids = []
    for i in range(k):
        z, w = nodes1(m[i], s[i], q=q, R=8.0); grids.append((z, w))
    Z = np.meshgrid(*[g[0] for g in grids], indexing="ij"); Wt = np.prod(np.meshgrid(*[g[1] for g in grids], indexing="ij"), axis=0)
    X = np.stack([z.ravel() - m[i] for i, z in enumerate(Z)], axis=1)
    Ci = np.linalg.inv(C); dens = np.exp(-0.5 * np.einsum("ni,ij,nj->n", X, Ci, X)) / np.sqrt((2 * np.pi) ** k * np.linalg.det(C))
    return np.sum(Wt.ravel() * dens * f(*[z.ravel() for z in Z]))


print("1. tokens vs 1-D quadrature")
for mm, ss in ((0.0, 1.3), (0.7, 0.9), (-1.1, 0.6)):
    t = {k: v.item() for k, v in tokens(np.array([mm]), np.array([ss * ss])).items()}
    G = E1(relu, mm, ss); G2 = E1(lambda z: relu(z) ** 2, mm, ss)
    M3 = E1(lambda z: relu(z) ** 3, mm, ss); M4 = E1(lambda z: relu(z) ** 4, mm, ss)
    k3 = M3 - 3 * G2 * G + 2 * G ** 3; k4 = M4 - 4 * M3 * G - 3 * G2 ** 2 + 12 * G2 * G ** 2 - 6 * G ** 4
    X = lambda z: (relu(z) - G) ** 2
    F1 = E1(lambda z: X(z) * (z - mm) / ss ** 2, mm, ss)                              # Stein: E dX = E X (z-m)/s^2
    F2 = E1(lambda z: X(z) * ((z - mm) ** 2 - ss ** 2) / ss ** 4, mm, ss)              # E d^2 X
    print(f"  m={mm:+.1f} s={ss}: dG {abs(t['G']-G):.1e} dG2 {abs(t['G2']-G2):.1e} dk3 {abs(t['k3']-k3):.1e} dk4 {abs(t['k4']-k4):.1e} dF1 {abs(t['F1']-F1):.1e} dF2 {abs(t['F2']-F2):.1e}")

print("2. pairs: residual of the two-token cov, coincident kappa3 and K^u formulas vs 2-D quadrature (expect ~c^3)")
ma, mb, sa, sb = 0.4, -0.3, 1.1, 0.8
ta = {k: v.item() for k, v in tokens(np.array([ma]), np.array([sa * sa])).items()}; tb = {k: v.item() for k, v in tokens(np.array([mb]), np.array([sb * sb])).items()}
for c in (0.2, 0.1, 0.05):
    C = np.array([[sa * sa, c], [c, sb * sb]]); mv = np.array([ma, mb])
    cov = Ek(lambda u, v: relu(u) * relu(v), mv, C) - ta["G"] * tb["G"]
    cov2 = ta["a1"] * tb["a1"] * c + 0.5 * ta["a2"] * tb["a2"] * c * c
    k3aab = Ek(lambda u, v: (relu(u) - ta["G"]) ** 2 * (relu(v) - tb["G"]), mv, C)
    k3f = ta["F1"] * tb["a1"] * c + 0.5 * ta["F2"] * tb["a2"] * c * c
    Xa = lambda u: (relu(u) - ta["G"]) ** 2; Xb = lambda v: (relu(v) - tb["G"]) ** 2
    Ku = Ek(lambda u, v: Xa(u) * Xb(v), mv, C) - E1(Xa, ma, sa) * E1(Xb, mb, sb) - 2 * cov ** 2
    Kf = ta["F1"] * tb["F1"] * c + (0.5 * ta["F2"] * tb["F2"] - 2 * ta["a1"] ** 2 * tb["a1"] ** 2) * c * c
    print(f"  c={c}: cov {cov:+.6f} resid {cov-cov2:+.1e} | k3_aab {k3aab:+.6f} resid {k3aab-k3f:+.1e} | K^u_ab {Ku:+.6f} resid {Ku-Kf:+.1e}")

print("3. distinct triple: cherry vs 3-D quadrature (expect residual ~C^3)")
m3 = np.array([0.3, -0.2, 0.5]); s3 = np.array([1.0, 0.9, 1.2]); t3 = [{k: v.item() for k, v in tokens(np.array([m3[i]]), np.array([s3[i] ** 2])).items()} for i in range(3)]
for c in (0.15, 0.075):
    C = np.diag(s3 ** 2); C[0, 1] = C[1, 0] = c; C[0, 2] = C[2, 0] = 0.8 * c; C[1, 2] = C[2, 1] = -0.6 * c
    k3 = Ek(lambda u, v, z: (relu(u) - t3[0]["G"]) * (relu(v) - t3[1]["G"]) * (relu(z) - t3[2]["G"]), m3, C, q=40)
    cherry = sum(t3[(h + 1) % 3]["a1"] * t3[(h + 2) % 3]["a1"] * t3[h]["a2"] * C[h, (h + 1) % 3] * C[h, (h + 2) % 3] for h in range(3))
    print(f"  c={c}: k3_abc {k3:+.4e} cherry {cherry:+.4e} resid {k3-cherry:+.1e}")

print("4. register readouts vs the represented tensor (sym = sum of the three placements)")
rng = np.random.default_rng(0); n = 6
T = rng.standard_normal((n, n)); TX = rng.standard_normal((n, n)); TD = rng.standard_normal((n, n)); a2 = rng.random(n)
Q = T * a2[None, :]; Ten = np.zeros((n, n, n))
for k in range(n):
    for P_, Q_ in ((TX, Q), (T, TD)):
        x, y = P_[:, k], Q_[:, k]
        Ten += np.einsum("a,b,c->abc", x, x, y) + np.einsum("a,b,c->abc", x, y, x) + np.einsum("a,b,c->abc", y, x, x)
k3, M = readouts([dict(T=T, TX=TX, TD=TD, a2=a2)])
print(f"  diag: {np.abs(k3 - np.einsum('aaa->a', Ten)).max():.1e}; (a,a,b): {np.abs(M - np.einsum('aab->ab', Ten)).max():.1e}")
