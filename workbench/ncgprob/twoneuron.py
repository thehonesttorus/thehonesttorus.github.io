# Independent numerical check of the two-neuron closed form of the pasted PDF: Z1, Z2 iid N(0,1), H = relu(Z), incoming
# pair source T_112 = T_121 = T_211 = t (score h_T = (t/2)(Z1^2 - 1) Z2), omitted-class response for the row w = (a, b):
#   Delta_T(a, b) = (3 t c / 2) a b [a^2 (1 + 4 c^2) + b^2 (1 - 4 c^2)],  c = 1/sqrt(2 pi),
# via the stated per-entry variations d cum(H1,H1,H1,H2) = t (3c/8 + 3c^3/2), d cum(H1,H2,H2,H2) = t (3c/8 - 3c^3/2),
# computed here by Gauss-Hermite quadrature of the score-weighted joint moments (exact to quadrature precision).
import numpy as np
from numpy.polynomial.hermite_e import hermegauss
x, wq = hermegauss(80); wq = wq / wq.sum()
Z1, Z2 = np.meshgrid(x, x, indexing="ij"); W2 = np.outer(wq, wq); H1, H2 = np.maximum(Z1, 0), np.maximum(Z2, 0); h = 0.5 * (Z1**2 - 1) * Z2   # t = 1
E = lambda f: np.sum(W2 * f); L = lambda f: np.sum(W2 * h * f)
m1, m2 = E(H1), E(H2); X1, X2 = H1 - m1, H2 - m2
# first variation of the joint cumulant cum(X1,X1,X1,X2) and cum(X1,X2,X2,X2) under the score tangent (means move too)
def dcum31(A, B):
    # cum(A,A,A,B) = E[A^3 B] - 3 E[A^2] E[A B]  for centred A, B; variation includes the recentering of A, B
    dmA, dmB = L(A), L(B)
    dE_A3B = L(A**3 * B) - 3 * dmA * E(A**2 * B) - dmB * E(A**3)
    dE_A2 = L(A**2) - 2 * dmA * E(A); dE_AB = L(A * B) - dmA * E(B) - dmB * E(A)
    return dE_A3B - 3 * (dE_A2 * E(A * B) + E(A**2) * dE_AB)
c = 1 / np.sqrt(2 * np.pi)
d1112 = dcum31(X1, X2); d1222 = dcum31(X2, X1)
print(f"d cum(H1,H1,H1,H2): quadrature {d1112:+.8f}, PDF {3*c/8 + 3*c**3/2:+.8f}")
print(f"d cum(H1,H2,H2,H2): quadrature {d1222:+.8f}, PDF {3*c/8 - 3*c**3/2:+.8f}")
for a, b in ((1.0, 1.0), (1.0, 0.5), (0.3, 1.0)):
    # omitted classes of kappa4(a X1 + b X2): (3+1) only here (no (2+1+1) or (1111) with two neurons): 4 a^3 b cum(1112) + 4 a b^3 cum(1222)
    num = 4 * a**3 * b * d1112 + 4 * a * b**3 * d1222; pdf = 1.5 * c * a * b * (a * a * (1 + 4 * c * c) + b * b * (1 - 4 * c * c))
    print(f"  w = ({a}, {b}): quadrature {num:+.8f}, PDF closed form {pdf:+.8f}")
