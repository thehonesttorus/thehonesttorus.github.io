# Exact checks behind the deck (Moreau) symmetry of the ReLU chain, by Gauss-Hermite quadrature at the Gaussian
# reference z ~ N(mu, s^2) (and a bivariate reference for the birth terms). No network.
#   (1) relu(z) - relu(-z) = z, so every derivative of order >= 2 of relu(z) and of its polar part relu(-z) agree:
#       the hinge coefficients c(1,k) = E[d^k relu(z)/dz^k], k >= 2, satisfy c(1,k)(alpha) = (-1)^k c(1,k)(-alpha)
#       (|c(1,k)| is even in alpha), while the transmission c(1,1) = Phi(alpha) is not (Phi(alpha) + Phi(-alpha) = 1).
#   (2) the cumulants of order >= 2 of an on-saturated unit's polar part v = relu(-z) equal those of the mirrored
#       off-saturated unit relu(z') with alpha' = -alpha, and decay with the Gaussian tail.
#   (3) a birth term whose nonlinearity sits at unit j: for (z_j, z_k) jointly Gaussian with y = relu(z),
#       kappa(y_j, y_j, y_k) = kappa(z_j, z_j, y_k) + 2 kappa(z_j, v_j, y_k) + kappa(v_j, v_j, y_k); the first term is
#       hub k's (the nonlinearity at k, j a linear leg), and the j-hub terms carry only v_j, so they vanish like the
#       tail of alpha_j on BOTH sides.
import numpy as np
from math import erf, sqrt, pi, exp

x, w = np.polynomial.hermite_e.hermegauss(200); w = w / w.sum()   # probabilists' Gauss-Hermite, E over N(0,1)
relu = lambda t: np.maximum(t, 0.0)
Phi = lambda a: 0.5 * (1 + erf(a / sqrt(2))); phi = lambda a: exp(-0.5 * a * a) / sqrt(2 * pi)


def hinge(alpha, s, k):
    # E[d^k relu / dz^k] for z ~ N(mu, s^2): k = 1 -> Phi; k >= 2 -> derivatives of the density at the hinge (Price)
    if k == 1:
        return Phi(alpha)
    # d^k relu = delta^(k-2): E[delta^(m)(z)] = (-1)^m p^(m)(0), p the N(mu, s^2) density; p^(m)(0) via Hermite
    m = k - 2; He = np.polynomial.hermite_e.HermiteE.basis(m)
    # p^(m)(z) = (-1)^m s^-(m+1) He_m((z - mu)/s) phi((z - mu)/s)  =>  E[delta^(m)] = s^-(m+1) He_m(-alpha) phi(alpha)
    return s ** (-(m + 1)) * He(-alpha) * phi(alpha)


print("(1) hinge coefficients c(1,k) under alpha -> -alpha (s = 0.7):")
for a in (0.3, 1.0, 2.5, 3.5):
    row = [f"k={k}: {hinge(a, 0.7, k):+.6e} vs (-1)^k x {(-1) ** k * hinge(-a, 0.7, k):+.6e}" for k in (2, 3, 4)]
    print(f"  alpha {a}: " + "  ".join(row) + f"   Phi(a)+Phi(-a) = {Phi(a) + Phi(-a):.12f}")
# check of c(1,2) and c(1,3) against finite differences in mu of the closed form E relu(z) = mu Phi(mu/s) + s phi(mu/s)
# (Price's theorem: d^k/dmu^k E relu(z) = E[d^k relu/dz^k])
for a in (0.5, 2.5):
    s = 0.7; mu = a * s; h = 1e-3
    E = lambda m_: m_ * Phi(m_ / s) + s * phi(m_ / s)
    d2 = (E(mu + h) - 2 * E(mu) + E(mu - h)) / h ** 2; d3 = (E(mu + 2 * h) - 2 * E(mu + h) + 2 * E(mu - h) - E(mu - 2 * h)) / (2 * h ** 3)
    print(f"  finite differences at alpha {a}: c(1,2) {d2:.6e} (exact {hinge(a, s, 2):.6e}), c(1,3) {d3:.4e} (exact {hinge(a, s, 3):.4e})")


def cums(vals):
    m = np.sum(w * vals); c = vals - m
    return [float(np.sum(w * c ** p)) for p in (2, 3, 4)]


print("(2) cumulants of the polar part relu(-z) at alpha vs relu(z') at -alpha (s = 1):")
for a in (1.0, 2.0, 2.5, 3.0):
    zp = a + x; zm = -a + x
    c_pol = cums(relu(-zp)); c_mir = cums(relu(zm))
    k4 = lambda c: c[2] - 3 * c[0] ** 2
    print(f"  alpha {a}: var {c_pol[0]:.3e}/{c_mir[0]:.3e}  k3 {c_pol[1]:.3e}/{c_mir[1]:.3e}  k4 {k4(c_pol):.3e}/{k4(c_mir):.3e}"
          f"   (phi(a)^2/phi(0)^2 = {phi(a) ** 2 / phi(0) ** 2:.1e})")

print("(3) birth term kappa(y_j, y_j, y_k) split at a bivariate Gaussian reference (rho = 0.4, s = 1):")
X1, X2 = np.meshgrid(x, x, indexing="ij"); W2 = np.outer(w, w)
rho = 0.4
for aj, ak in ((3.0, 0.2), (2.5, -0.4), (-2.5, 0.2), (0.3, 0.2)):
    zj = aj + X1; zk = ak + rho * X1 + sqrt(1 - rho * rho) * X2
    yj, yk, vj = relu(zj), relu(zk), relu(-zj)
    E = lambda f: float(np.sum(W2 * f)); cen = lambda f: f - E(f)
    k3 = lambda a_, b_, c_: E(cen(a_) * cen(b_) * cen(c_))
    tot = k3(yj, yj, yk); hub_k = k3(zj, zj, yk); hub_j = 2 * k3(zj, vj, yk) + k3(vj, vj, yk)
    print(f"  alpha_j {aj:+.1f}, alpha_k {ak:+.1f}: total {tot:+.4e} = hub k (z_j linear) {hub_k:+.4e} + hub j (polar) {hub_j:+.4e}"
          f"  [check {tot - hub_k - hub_j:+.1e}]")
