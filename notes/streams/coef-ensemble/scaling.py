"""Width scaling of the closure ladder and of the ensemble coefficients (reads coef_table.py's summary.json).

    python scaling.py [results/summary.json] > results/scaling.txt

eps: per layer and model, least squares of log eps on log n over the widths present (noise-corrected cross-evaluated
eps on the pair MLPs, and the within-atlas mean over all MLPs), giving eps ~ A n^-p and the extrapolation to n = 1024.
Coefficients: per layer and basis term, c(n) = c_inf + a / sqrt(n) (the first-order corrections are O(rho) ~ n^-1/2),
fitted over the widths, giving the table value at n = 1024.
"""
import json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "results", "summary.json")
S = json.load(open(path))
widths = sorted(int(w) for w in S)
MODELS = ["wick", "leg", "legR", "own", "ensD", "ens", "ensR", "ensRD"]
NAMES = ["B0 k3z", "B1 D21+C", "B2 w3D21", "B3 D3", "B4 rho3", "B5 K22", "B6 K211"]
N_EXT = 1024


def powfit(ns, es):
    ns, es = np.array(ns, float), np.array(es, float)
    ok = es > 0
    if ok.sum() < 2:
        return np.nan, np.nan
    p, c = np.polyfit(np.log(ns[ok]), np.log(es[ok]), 1)
    return -p, float(np.exp(c) * N_EXT ** p)


L1 = len(S[str(widths[0])]["noise"])
print(f"widths {widths}; MLPs per width " + ", ".join(f"{w}: {S[str(w)]['n_mlps']} ({S[str(w)]['n_pairs']} pairs)" for w in widths))
def floorfit(ns, es):
    """eps = a + b / sqrt(n): a floor plus the O(rho) term; the conservative extrapolation"""
    ns, es = np.array(ns, float), np.array(es, float)
    ok = np.isfinite(es)
    if ok.sum() < 2:
        return np.nan
    X = np.stack([np.ones(ok.sum()), 1 / np.sqrt(ns[ok])], 1)
    (a, b), *_ = np.linalg.lstsq(X, es[ok], rcond=None)
    return float(a + b / np.sqrt(N_EXT))


for key, label in [("eps_rep", "representation error eps_rep (pair MLPs; target and model-input MC noise subtracted)"),
                   ("eps_raw", "within-atlas eps (all MLPs, includes MC noise)")]:
    print(f"\n## {label}: per width, then fit eps ~ n^-p and the value at n = {N_EXT}; last column the floor model a + b/sqrt(n) at {N_EXT}")
    for m in MODELS:
        print(f"\n### {m}")
        print(" l | " + " ".join(f"n={w:<4d}" for w in widths) + " |     p  eps(1024) floor(1024)")
        for l in range(L1):
            es = [S[str(w)][key][m][l] for w in widths]
            p, e = powfit(widths, es)
            print(f"{l:>2} | " + " ".join(f"{v:6.3f}" for v in es) + f" | {p:5.2f}  {e:7.4f}  {floorfit(widths, es):7.4f}")
        for lo, hi in [(1, 3), (4, 9), (10, 14)]:
            es = [np.sqrt(np.mean(np.square(S[str(w)][key][m][lo:hi + 1]))) for w in widths]
            p, e = powfit(widths, es)
            print(f"rms {lo}-{hi} | " + " ".join(f"{v:6.3f}" for v in es) + f" | {p:5.2f}  {e:7.4f}  {floorfit(widths, es):7.4f}")

# same-seed comparison: only the MLP seeds that have a pair atlas at every width (removes MLP-sampling noise between
# widths, at the cost of fewer MLPs)
common = sorted(set.intersection(*[set(S[str(w)]["pair_seeds"]) for w in widths]))
print(f"\n## eps_rep on the MLP seeds with a pair at every width: {common}; band rms of the per-seed mean, fit n^-p, value at {N_EXT}")
for m in MODELS:
    rows = []
    for w in widths:
        sd = S[str(w)]["seeds"]
        E = np.array(S[str(w)]["eps_rep_all"][m], float)
        rows.append(np.nanmean(E[[sd.index(c) for c in common]], 0))
    for lo, hi in [(1, 3), (4, 9), (10, 14)]:
        es = [np.sqrt(np.mean(r[lo:hi + 1] ** 2)) for r in rows]
        p, e = powfit(widths, es)
        print(f"{m:>6} {lo}-{hi} | " + " ".join(f"{v:6.3f}" for v in es) + f" | {p:5.2f}  {e:7.4f}  {floorfit(widths, es):7.4f}")

for tab in ["tensor", "D21", "tensor_reg", "D21_reg"]:
    print(f"\n## ensemble coefficients ({tab} fit) per width, and c(n) = c_inf + a/sqrt(n) evaluated at n = {N_EXT}")
    nm = NAMES[:6] + (["B6 K211"] if "reg" not in tab else ["B6 uC"])
    for k in range(7):
        print(f"\n### {nm[k]}")
        print(" l | " + " ".join(f"n={w:<4d}" for w in widths) + " | c(1024)  c_inf  (across-MLP sd at each width)")
        for l in range(L1):
            cs = [S[str(w)]["coef"][tab][l][k] for w in widths]
            sds = [S[str(w)]["coef_sd"][tab][l][k] for w in widths]
            X = np.stack([np.ones(len(widths)), 1 / np.sqrt(widths)], 1)
            (cinf, a), *_ = np.linalg.lstsq(X, np.array(cs), rcond=None)
            print(f"{l:>2} | " + " ".join(f"{v:+6.2f}" for v in cs) + f" | {cinf + a / np.sqrt(N_EXT):+6.2f} {cinf:+6.2f}  ("
                  + " ".join(f"{v:4.2f}" for v in sds) + ")")

# the deliverable data file: ensemble tables per width and the n = 1024 extrapolation (tensor-space fit, true and
# regenerated (2,1,1) slice), layers 0..L-2 (row l = closure of kappa3(a_l) feeding D21(l+1)); layer-0 B0..B3, B5, B6
# multiply population-zero inputs (z_0 is exactly Gaussian) and are set to their leg-partition values
LEG = [1.0, 3.0, 3.0, 1.0, 1.0, 1.5, 1.5]
tables = {"basis": NAMES, "leg_partition": LEG, "per_width": {}, "n1024": {}}
for tab in ["tensor", "tensor_reg"]:
    C = np.array([S[str(w)]["coef"][tab] for w in widths])          # (widths, layers, 7)
    X = np.stack([np.ones(len(widths)), 1 / np.sqrt(widths)], 1)
    coef, *_ = np.linalg.lstsq(X, C.reshape(len(widths), -1), rcond=None)
    ext = (coef[0] + coef[1] / np.sqrt(N_EXT)).reshape(C.shape[1:])
    for arr in [ext] + [C[i] for i in range(len(widths))]:
        arr[0] = [LEG[k] if k != 4 else arr[0][4] for k in range(7)]
    tables["n1024"][tab] = np.round(ext, 3).tolist()
    tables["per_width"][tab] = {str(w): np.round(C[i], 3).tolist() for i, w in enumerate(widths)}
json.dump(tables, open(os.path.join(HERE, "results", "coef_table_n1024.json"), "w"), indent=1)
