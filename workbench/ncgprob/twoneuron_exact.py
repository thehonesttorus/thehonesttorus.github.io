# Exact check of the PDF's two-neuron witness with 1-D adaptive quadrature (the two coordinates are independent, so every
# expectation factorises into half-normal integrals; the quadrature is on smooth pieces split at the ReLU kink).
import numpy as np
from scipy.integrate import quad
ph = lambda z: np.exp(-z * z / 2) / np.sqrt(2 * np.pi)
def E(f):   # E f(Z), Z ~ N(0,1), split at 0
    return quad(lambda z: f(z) * ph(z), -np.inf, 0, epsabs=1e-13, epsrel=1e-13)[0] + quad(lambda z: f(z) * ph(z), 0, np.inf, epsabs=1e-13, epsrel=1e-13)[0]
c = 1 / np.sqrt(2 * np.pi); m = E(lambda z: max(z, 0)); U = lambda z: max(z, 0) - m
v = E(lambda z: U(z)**2); k3 = E(lambda z: U(z)**3)
# source score h = (t/2)(Z1^2 - 1) Z2 with t = 1; first variation of cum(X1,X1,X1,X2) with recentering, as in twoneuron.py, factorised
s1 = lambda f: E(lambda z: 0.5 * (z * z - 1) * f(z))      # the Z1 factor of the score against a function of Z1
s2 = lambda f: E(lambda z: z * f(z))                      # the Z2 factor
# L[g(Z1) k(Z2)] = s1(g) s2(k);  E[g(Z1) k(Z2)] = E(g) E(k)
dm1 = s1(lambda z: U(z)) * s2(lambda z: 1.0); dm2 = s1(lambda z: 1.0) * s2(lambda z: U(z))
L = lambda g, k: s1(g) * s2(k); EE = lambda g, k: E(g) * E(k)
U1, U2 = U, U
def dcum31():
    dE_A3B = L(lambda z: U1(z)**3, U2) - 3 * dm1 * EE(lambda z: U1(z)**2, U2) - dm2 * EE(lambda z: U1(z)**3, lambda z: 1.0)
    dE_A2 = L(lambda z: U1(z)**2, lambda z: 1.0) - 2 * dm1 * EE(U1, lambda z: 1.0)
    dE_AB = L(U1, U2) - dm1 * EE(lambda z: 1.0, U2) - dm2 * EE(U1, lambda z: 1.0)
    return dE_A3B - 3 * (dE_A2 * EE(U1, U2) + EE(lambda z: U1(z)**2, lambda z: 1.0) * dE_AB)
def dcum13():
    dE_AB3 = L(U1, lambda z: U2(z)**3) - dm1 * EE(lambda z: 1.0, lambda z: U2(z)**3) - 3 * dm2 * EE(U1, lambda z: U2(z)**2)
    dE_B2 = L(lambda z: 1.0, lambda z: U2(z)**2) - 2 * dm2 * EE(lambda z: 1.0, U2)
    dE_AB = L(U1, U2) - dm1 * EE(lambda z: 1.0, U2) - dm2 * EE(U1, lambda z: 1.0)
    return dE_AB3 - 3 * (dE_B2 * EE(U1, U2) + EE(lambda z: 1.0, lambda z: U2(z)**2) * dE_AB)
d1112, d1222 = dcum31(), dcum13()
print(f"d cum(H1,H1,H1,H2): exact quad {d1112:+.10f}, PDF 3c/8 + 3c^3/2 = {3*c/8 + 3*c**3/2:+.10f}")
print(f"d cum(H1,H2,H2,H2): exact quad {d1222:+.10f}, PDF 3c/8 - 3c^3/2 = {3*c/8 - 3*c**3/2:+.10f}")
for a, b in ((1.0, 1.0), (1.0, 0.5), (0.3, 1.0)):
    num = 4 * a**3 * b * d1112 + 4 * a * b**3 * d1222; pdf = 1.5 * c * a * b * (a * a * (1 + 4 * c * c) + b * b * (1 - 4 * c * c))
    print(f"  w = ({a}, {b}): exact {num:+.10f}, PDF {pdf:+.10f}, rel diff {abs(num-pdf)/abs(pdf):.1e}")
