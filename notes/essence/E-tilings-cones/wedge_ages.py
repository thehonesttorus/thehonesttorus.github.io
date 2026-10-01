"""Annealed wedge model of the age profile of old third-order content in D21 (compare region REPORT N6, faces T3).
Energy of the D21 contribution of source layer s at target t (FC's Wick term, incoherent over atoms r):
  E_s(t) ~ B_s * prod_{u=s+1}^{t-1} lambda_u^3,   lambda_u = 1/2 + arcsin(rho_u)/pi  (one wedge factor per leg),
  B_s = sigma_s^8 E[w2^2 Phi^4],  w2 = phi(t)/sigma (facet density), t ~ N(0, rho/(1-rho)), sigma^2 = 2(1-rho).
Ages are mutually orthogonal (F8.3), so dropping ages > k leaves relative D21 error sqrt(1 - share(ages <= k))."""
import numpy as np
from math import erf
L = 16
f = lambda r: (np.sqrt(1 - r * r) + (np.pi - np.arccos(r)) * r) / np.pi
rho = [0.0]
for _ in range(L - 1):
    rho.append(f(rho[-1]))
rho = np.array(rho); lam = 0.5 + np.arcsin(rho) / np.pi
tau2 = rho / (1 - rho); sig2 = 2 * (1 - rho)
x, w = np.polynomial.hermite_e.hermegauss(80); w /= w.sum()
Phi = lambda v: 0.5 * (1 + np.vectorize(erf)(v / np.sqrt(2)))
phi = lambda v: np.exp(-v * v / 2) / np.sqrt(2 * np.pi)
B = np.array([sig2[s] ** 4 * np.sum(w * (phi(np.sqrt(tau2[s]) * x) ** 2 / sig2[s]) * Phi(np.sqrt(tau2[s]) * x) ** 4)
              for s in range(L)])
print('birth strength B_s/B_0 :', np.round(B / B[0], 4))
print('three-leg transport lambda^3:', np.round(lam ** 3, 3))
for tgt in [8, 11, 14]:
    E = np.array([B[s] * np.prod(lam[s + 1:tgt] ** 3) for s in range(tgt)]); E /= E.sum()
    age = tgt - np.arange(tgt)
    print(f'target {tgt}: age shares (age 1..{tgt})', np.round(E[::-1], 3))
    print(f'   keep ages<=2: D21 rel err {np.sqrt(1 - E[age <= 2].sum()):.2f}; keep ages<=4: {np.sqrt(1 - E[age <= 4].sum()):.2f}')
print('region N6 measured at layer 14 (MLP 0, inside FC): ages<=2 -> 0.89, ages<=4 -> 0.75')
