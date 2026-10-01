"""Markdown tables for REPORT.md from results/summary.json (coef_table.py analyse)."""
import json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
S = json.load(open(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "results", "summary.json")))
widths = sorted(int(w) for w in S)
NAMES = ["B0 κ3(z)", "B1 D21⊗C", "B2 w3 D21", "B3 D3⊗C⊗C", "B4 ρ³", "B5 K22⊗C", "B6 κ4(2,1,1)"]
LEG = [1, 3, 3, 1, 1, 1.5, 1.5]
BANDS = [(1, 3), (4, 9), (10, 14)]

print("### Held-out representation error ε_rep of D21(l+1) (noise-corrected, mean over pair MLPs; per layer)\n")
models = ["wick", "leg", "legR", "own", "ownD", "ens", "ensD", "ensR", "ensRD"]
for w in widths:
    s = S[str(w)]
    print(f"\nwidth {w} ({s['n_mlps']} MLPs, {s['n_pairs']} pairs)\n")
    print("| l | " + " | ".join(models) + " |")
    print("|---" * (len(models) + 1) + "|")
    for l in range(1, len(s["noise"])):
        print(f"| {l} | " + " | ".join(f"{s['eps_rep'][m][l]:.3f}" for m in models) + " |")

print("\n### Band summary (rms over layers of the mean ε_rep)\n")
print("| width | " + " | ".join(f"{m} {a}-{b}" for m in ("leg", "own", "ens", "ensR") for a, b in BANDS) + " |")
print("|---" * (1 + 4 * len(BANDS)) + "|")
for w in widths:
    s = S[str(w)]
    cells = []
    for m in ("leg", "own", "ens", "ensR"):
        e = np.array(s["eps_rep"][m])
        cells += [f"{np.sqrt(np.nanmean(e[a:b + 1] ** 2)):.3f}" for a, b in BANDS]
    print(f"| {w} | " + " | ".join(cells) + " |")

print("\n### Ensemble coefficients (tensor-space fit over all MLPs of the width) ± across-MLP sd of per-MLP fits\n")
for k in range(7):
    print(f"\n{NAMES[k]} (leg-partition value {LEG[k]})\n")
    print("| l | " + " | ".join(f"n = {w}" for w in widths) + " |")
    print("|---" * (len(widths) + 1) + "|")
    for l in range(1, len(S[str(widths[0])]["noise"])):
        print(f"| {l} | " + " | ".join(f"{S[str(w)]['coef']['tensor'][l][k]:+.2f} ± {S[str(w)]['coef_sd']['tensor'][l][k]:.2f}" for w in widths) + " |")
