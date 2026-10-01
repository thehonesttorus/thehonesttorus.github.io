"""Build the §4.6 table from oracle2nd rows and the computed variants on the same w128 networks."""
import json, glob, collections, numpy as np, os
HERE = os.path.dirname(os.path.abspath(__file__))
rows = [json.loads(l) for f in glob.glob(os.path.join(HERE, "results_live", "oracle2nd_*.jsonl")) for l in open(f)]
orc = collections.defaultdict(dict)
for r in rows: orc[(r["mode"], r["r"])][r["mlp"]] = r["raw"]
comp = collections.defaultdict(dict)
for f in ("w128_d16.jsonl", "w128_d16_nc.jsonl"):
    for r in map(json.loads, open(os.path.join(HERE, "results_live", f))): comp[r["variant"]][r["mlp"]] = r["raw"]
mlps = sorted(orc[("full2", None)])
fmt = lambda d: " / ".join(f"{d[i]:.2e}" for i in mlps) + f" | {np.exp(np.mean(np.log([d[i] for i in mlps]))):.2e}"
print(f"\n\n| (w128, networks {', '.join(map(str, mlps))}) | raw per network | geometric mean |\n|---|---|---|")
for lab, d in [("Gaussian closure", comp["gauss"]), ("computed, exact first order (all pairs)", comp["full"]),
               ("computed, best co-state A3gsl_nc", comp["A3gsl_nc"]),
               ("oracle: true κ₃, first-order injection", orc[("first", None)]),
               ("oracle: true κ₃ + joint κ₄ (second order)", orc[("full2", None)]),
               ("oracle: κ₄ slices rank 1 off-diagonal + exact diagonal", orc[("full2", 1)]),
               ("oracle: κ₄ slices rank 8 + diagonal", orc[("full2", 8)]),
               ("oracle: κ₄ slices rank 32 + diagonal", orc[("full2", 32)]),
               ("oracle: κ₄ diagonal only", orc[("k4diag", None)])]:
    print(f"| {lab} | {fmt(d)} |")
