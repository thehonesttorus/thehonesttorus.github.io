"""Contact decomposition of the layer-two third cumulant (notes/stage10, Theorem: cumulants are walk sums).
  python scripts/verify_tokens_contacts.py"""
import numpy as np
from math import sqrt, pi, factorial
from scipy.integrate import quad
rng = np.random.default_rng(11)
phi0 = 1 / sqrt(2 * pi)
# Hermite coefficients E[f He_k] on the half line for f = relu and F = relu^2 - 2 phi0 relu
He = np.polynomial.hermite_e.HermiteE
def coef(fun, K):
    return np.array([quad(lambda s: fun(s) * He.basis(k)(s) * np.exp(-s * s / 2) / sqrt(2 * pi), 0, 40, limit=400)[0] for k in range(K)])
K = 60
c = coef(lambda s: s, K)                     # relu
a = coef(lambda s: s * s, K) - 2 * phi0 * c  # relu^2 - 2 phi0 relu
k2 = lambda r: sum(a[k] * c[k] * r ** k / factorial(k) for k in range(1, K))     # kappa_3(h,h,h') per unit norms
print("Hermite coefficients of relu c_0..c_6:", np.round(c[:7], 4))
for n in (24, 96):
    W1 = rng.standard_normal((n, n)) * sqrt(2 / n); w = rng.standard_normal(n) * sqrt(2 / n)
    nw = np.linalg.norm(W1, axis=1); R = (W1 @ W1.T) / np.outer(nw, nw); om = w * nw   # om_a = W2(c,a)|w_a|
    k3_self = np.sum(om ** 3) * (0.5 * phi0 + 2 * phi0 ** 3)
    off = ~np.eye(n, dtype=bool)
    k3_pair = 3 * np.sum((om ** 2)[:, None] * om[None, :] * np.vectorize(k2)(np.where(off, R, 0)) * off)
    # all-distinct, leading diagrams: hub paths (degree 2 at the hub, gates at the arms) + triangles
    Ro = R * off; T1 = om @ Ro                                              # sum_b om_b R_ab
    paths = 3 * c[1] ** 2 * c[2] * np.sum(om * (T1 ** 2 - (om[:, None] ** 0 * (om ** 2)[None, :] * Ro ** 2).sum(1)))
    tri = c[2] ** 3 * np.einsum('a,b,e,ab,be,ea->', om, om, om, Ro, Ro, Ro)
    N = 1 << 23; mom = np.zeros(3)
    for _ in range(16):
        x = rng.standard_normal((N // 16, n)); z = np.maximum(x @ W1.T, 0) @ w; mom += [z.sum(), (z ** 2).sum(), (z ** 3).sum()]
    mom /= N; mc = mom[2] - 3 * mom[1] * mom[0] + 2 * mom[0] ** 3
    distinct_mc = mc - k3_self - k3_pair
    print(f"n={n}: kappa_3 MC {mc:.5f} = self (one neuron) {k3_self:.5f} + pairs (two neurons, one repeated) {k3_pair:.5f} + distinct {distinct_mc:.5f}")
    print(f"       distinct part, leading diagrams: hub paths {paths:.5f} + triangles {tri:.5f} = {paths + tri:.5f}")
