"""Paired comparison of lab_run variants against a baseline on the same MLPs.

  python compare.py --res RESDIR --base TAG [TAG ...]

Per variant: mean final MSE, mean (MSE - truth noise) = 'raw', the paired difference vs the baseline
(mean over MLPs, its s.e. across MLPs, relative change of raw), C/B, max residual, and the truth-noise
part of the s.e. of the paired difference (2 sqrt(sum_i dp_i^2 noise) / n per MLP, combined).
"""
import argparse, json, os
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--res", required=True)
ap.add_argument("--base", required=True)
ap.add_argument("tags", nargs="*")
ap.add_argument("--md", action="store_true")
a = ap.parse_args()


def load(tag):
    j = json.load(open(os.path.join(a.res, f"{tag}.json")))
    p = np.load(os.path.join(a.res, "preds", f"{tag}.npz"))
    return j, p["preds"], list(p["mlps"])


jb, pb, mb = load(a.base)
rb = {r["mlp"]: r for r in jb["per_mlp"]}
base_raw = np.mean([rb[m]["adj"] for m in mb])
hdr = "| variant | MLPs | final MSE | raw (MSE - noise) | Δ vs base (paired) | s.e. | Δ raw % | C/B | max residual s |"
print(hdr); print("|" + "---|" * 9)
for tag in [a.base] + a.tags:
    try:
        j, p, mm = load(tag)
    except FileNotFoundError:
        print(f"| {tag} | missing |"); continue
    r = {x["mlp"]: x for x in j["per_mlp"]}
    common = [m for m in mm if m in rb and r[m].get("ok") and rb[m].get("ok")]
    d = np.array([r[m]["mse"] - rb[m]["mse"] for m in common])
    raw = np.mean([r[m]["adj"] for m in common])
    braw = np.mean([rb[m]["adj"] for m in common])
    se = d.std(ddof=1) / np.sqrt(len(d)) if len(d) > 1 else float("nan")
    print(f"| {tag} | {len(common)} | {np.mean([r[m]['mse'] for m in common]):.4e} | {raw:.4e} | {d.mean():+.3e} | {se:.1e} | "
          f"{100 * d.mean() / braw:+.1f} | {j['cb']:.4f} | {j['residual_max']:.3f} |")
