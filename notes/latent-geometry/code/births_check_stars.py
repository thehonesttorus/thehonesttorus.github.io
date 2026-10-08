# Star check: the odd part of kappa_4 in c3 is star(c3) + O(c3^3); extract the linear coefficient from c3 and 2 c3
# (common random numbers): linear = (8 odd(c3) - odd(2 c3)) / 6.
import numpy as np
exec(open("births_check.py").read().split("c1 = rng.uniform")[0])
c1 = rng.uniform(0.3, 1.0, n); c2 = rng.uniform(-0.8, 0.8, n); c3 = rng.uniform(-0.3, 0.3, n); z2 = np.zeros(n)
stp, _, sgp, _, _, _ = trees(c1, z2, c3)
odd = {}
for f in (1, 2):
    kp, Kp = mc(c1, z2, f * c3, seed=7); km, Km = mc(c1, z2, -f * c3, seed=7); odd[f] = ((kp - km) / 2, (Kp - Km) / 2)
lin_k = (8 * odd[1][0] - odd[2][0]) / 6; lin_K = (8 * odd[1][1] - odd[2][1]) / 6
off = ~np.eye(m, dtype=bool)
print("stars, kappa_4 diagonal: MC linear part", np.round(lin_k, 4), " star", np.round(sgp, 4))
print("stars, (3,1) offdiag: MC linear part", np.round(lin_K[off], 4), "\n    star", np.round(stp[off], 4))
