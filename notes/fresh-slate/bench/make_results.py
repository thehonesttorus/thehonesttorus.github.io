"""Regenerate RESULTS.md: the two calibration baselines on every set in sets/ (python make_results.py)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bench, eval_q  # noqa: E402

lines = ["# Baselines on the bench sets (calibration only, not designs)", "",
         "gauss = Gaussian covariance closure with linearised cross-covariance (`eval_q.baseline_gauss`, float64);",
         "mc = plain Monte Carlo with 6554 samples (the 0.1-floor budget at n = 1024: 0.1 B / (2 L n²)).",
         "raw = final-layer MSE − truth noise; mean over the set's MLPs (± s.e. across MLPs).", "",
         "| set | MLPs | N | truth noise | gauss raw | mc raw | gauss all-layer |", "|---|---|---|---|---|---|---|"]
res = {"gauss": [], "mc": []}
for s in bench.list_sets():
    g = eval_q.eval_set(eval_q.baseline_gauss, s["name"], verbose=False)
    m = eval_q.eval_set(eval_q.baseline_mc, s["name"], verbose=False)
    res["gauss"].append(g); res["mc"].append(m)
    lines.append(f"| {s['name']} | {s['n_mlps']} | {s['N']:.0e} | {g['noise']:.1e} | {g['raw']:.3e} ± {g['raw_se']:.1e} | "
                 f"{m['raw']:.3e} ± {m['raw_se']:.1e} | {g['all_layers']:.3e} |")
f = eval_q.scaling_fit(res["gauss"], units_1024=32.0)
if "p" in f:
    lines += ["", f"gauss width fit (depth-16 sets): raw ∝ n^-{f['p']:.2f}, extrapolated raw(1024) = {f['raw_1024']:.2e}"
              + (f"; measured at 1024: {f['raw_1024_measured']:.2e}" if "raw_1024_measured" in f else "")]
open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "RESULTS.md"), "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
