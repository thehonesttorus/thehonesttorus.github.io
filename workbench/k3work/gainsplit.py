# Split the chain's kappa3 / kappa4 errors into the gain-mode template and the remainder (ray-compiler note section 6).
#   python gainsplit.py NET MCPREFIX DUMP
# Per layer, against Monte Carlo truth (cross-half for the noise floor): the template amplitudes t (least squares on
# 6 mu var for the kappa3 diagonal, 12 var^2 for the kappa4 diagonal), the truth's non-gain remainder |X - t T| / |X|,
# and the chain's error e = chain - truth split into its component along the template and the rest, readout-weighted
# by the Edgeworth coefficients c3, c4 of the mean (what the means actually read).
import sys, math, numpy as np
net, pre, dp = int(sys.argv[1]), sys.argv[2], sys.argv[3]
T = {h: np.load(f"{pre}_{h}.npz") for h in ("full", "h0", "h1")}; D = np.load(dp)
L = T["full"]["mu"].shape[0]
phi = lambda a: np.exp(-a * a / 2) / math.sqrt(2 * math.pi)
print(f"net {net}: kappa3 diag (template 6 mu var) and kappa4 diag (12 var^2); readout-weighted (c3, c4) energies")
print("  layer | truth: non-gain fraction k3 k4 | chain error: |gain part| |rest| (k3) | |gain part| |rest| (k4) | noise floor k3 k4")
for l in range(2, L):
    mu = T["full"]["mu"][l].astype(np.float64); var = T["full"]["var"][l].astype(np.float64); s = np.sqrt(var); a = mu / s
    c3 = (-1) ** 0 * (1 / s) ** 2 * (-a) * phi(a) / 6.0 * 0 + (np.polynomial.hermite_e.HermiteE.basis(1)(-a) * phi(a) / s ** 2) / 6.0
    c4 = (np.polynomial.hermite_e.HermiteE.basis(2)(-a) * phi(a) / s ** 3) / 24.0
    out = []
    for key, tmpl, c in (("k3", 6 * mu * var, c3), ("k4", 12 * var * var, c4)):
        tru = T["full"][key][l].astype(np.float64)
        dk = {"k3": "D3", "k4": "g4row"}[key]
        if f"{dk}_{l}" not in D.files:
            out.append(None); continue
        ch = D[f"{dk}_{l}"].astype(np.float64)
        tt = float(tru @ tmpl / (tmpl @ tmpl)); ng = np.linalg.norm(tru - tt * tmpl) / np.linalg.norm(tru)
        e = ch - tru; we = c * e; wt = c * tmpl
        ge = float(we @ wt / (wt @ wt)) * wt; rest = we - ge
        h0, h1 = T["h0"][key][l].astype(np.float64), T["h1"][key][l].astype(np.float64)
        nf = np.linalg.norm(c * (h0 - h1)) / 2
        out.append((ng, np.linalg.norm(ge), np.linalg.norm(rest), nf))
    if None in out:
        continue
    (g3, a3, r3, n3), (g4, a4, r4, n4) = out
    print(f"   {l:2d}  | {g3:.3f} {g4:.3f}                    | {a3:.2e} {r3:.2e}         | {a4:.2e} {r4:.2e}          | {n3:.1e} {n4:.1e}")
