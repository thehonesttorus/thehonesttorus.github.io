# Effective dimension of the hidden state's fluctuation by layer (ray-compiler note section 6).
#   python effdim.py NET MCPREFIX
# Spectrum of the pre-activation covariance C^z_l and the post-activation covariance C^y_l: participation ratio
# (tr C)^2 / tr C^2, the number of directions holding 90 / 99 / 99.9% of the trace, and the power-law decay exponent of
# the eigenvalues over ranks 10-300 (lambda_k ~ k^-p).
import sys, numpy as np
net, pre = int(sys.argv[1]), sys.argv[2]
F = np.load(f"{pre}_full.npz")
L = F["mu"].shape[0]
print(f"net {net}: covariance spectra by layer (pre-activation z | post-activation y)")
print("  layer | z: PR   k90  k99  k999  decay p | y: PR   k90  k99  k999  decay p | y: mean^2/var (rank-one mean share)")
for l in range(L):
    row = []
    for key in ("cov", "cov_y"):
        C = F[key][l].astype(np.float64); C = 0.5 * (C + C.T)
        ev = np.clip(np.linalg.eigvalsh(C)[::-1], 0, None); tr = ev.sum(); cs = np.cumsum(ev) / tr
        pr = tr * tr / float(ev @ ev)
        ks = [int(np.searchsorted(cs, q) + 1) for q in (0.9, 0.99, 0.999)]
        k = np.arange(10, 300); p = -np.polyfit(np.log(k), np.log(ev[9:299] + 1e-300), 1)[0]
        row.append((pr, *ks, p))
    my = F["mu_y"][l].astype(np.float64); vy = F["var_y"][l].astype(np.float64)
    (pz, a1, a2, a3, dz), (py, b1, b2, b3, dy) = row
    print(f"   {l:2d}  | {pz:6.1f} {a1:4d} {a2:4d} {a3:4d}  {dz:5.2f}  | {py:6.1f} {b1:4d} {b2:4d} {b3:4d}  {dy:5.2f}  | {np.mean(my ** 2) / np.mean(vy):.2f}")
