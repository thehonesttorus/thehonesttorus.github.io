"""Error law MSE vs eps(D21) at width 128.

(1) injected: results/eps2/N<eps>s<seed>@<base>_mlp<i>.json  vs the unperturbed base  -> dMSE(final) = k eps^2
    (eps = relative rms of a random perturbation of D21(z) added at EVERY layer, as in 504aldo's F71f law)
(2) teacher forcing: F<l>@base (kappa3(a_l) := atlas), FD<l>@base (only its all-distinct part)
(3) across variants: results/eps/*.json with per-layer eps of D21 vs the atlas (noise-corrected if a pair-atlas noise
    file results/noise_mlp<i>.json exists) against final MSE
    python law.py [--base fit:atlas]
"""
import argparse, glob, json, os, re
import numpy as np


def load(d):
    out = {}
    for f in glob.glob(os.path.join(d, "*.json")):
        j = json.load(open(f)); out[(j["variant"], j["mlp"])] = j
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--base", default="fit:atlas")
    a = ap.parse_args()
    R = load("results/eps2")
    mlps = sorted({m for (_, m) in R})
    print(f"## (1) injected D21 noise at every layer, base {a.base}\n")
    print("| eps | " + " | ".join(f"dMSE final mlp{m}" for m in mlps) + " | k = dMSE/eps^2 (mean) | dMSE mean-over-layers / eps^2 |")
    print("|---|" + "---|" * len(mlps) + "---|---|")
    ks_all = []
    for eps in sorted({float(re.match(r"N([0-9.]+)s", v).group(1)) for (v, _) in R if v.startswith("N")}):
        d, dl = [], []
        for m in mlps:
            b = R.get((a.base, m))
            if b is None:
                continue
            reps = [R[(v, mm)] for (v, mm) in R if mm == m and v.startswith(f"N{eps}s") and v.endswith("@" + a.base)]
            if not reps:
                continue
            d.append(np.mean([r["final_mse"] - b["final_mse"] for r in reps]))
            dl.append(np.mean([np.mean(r["per_layer_mse"]) - np.mean(b["per_layer_mse"]) for r in reps]))
        k = np.mean(d) / eps ** 2
        ks_all.append((eps, k))
        print(f"| {eps} | " + " | ".join(f"{x:.2e}" for x in d) + f" | {k:.2e} | {np.mean(dl)/eps**2:.2e} |")
    print("\n## (2) teacher forcing (kappa3(a_l) replaced by the atlas's)\n")
    print("| variant | " + " | ".join(f"final mlp{m}" for m in mlps) + " | ratio to base (mean) |")
    print("|---|" + "---|" * len(mlps) + "---|")
    for v in sorted({v for (v, _) in R if v.startswith("F")}, key=lambda s: (s[:2] == "FD", int(re.search(r"\d+", s).group()))):
        vals = [R[(v, m)]["final_mse"] for m in mlps if (v, m) in R]
        base = [R[(a.base, m)]["final_mse"] for m in mlps if (v, m) in R]
        print(f"| {v} | " + " | ".join(f"{x:.2e}" for x in vals) + f" | {np.mean(np.array(vals)/np.array(base)):.2f} |")
    print(f"| {a.base} (base) | " + " | ".join(f"{R[(a.base, m)]['final_mse']:.2e}" for m in mlps if (a.base, m) in R) + " | 1 |")
    E = load("results/eps")
    print("\n## (3) across variants (teacher-forced kappa4): mean eps(D21) over layers 1-15 vs final MSE\n")
    print("| variant | mlp | mean eps D21 (raw) | mean eps D21 (noise-corr.) | final MSE |")
    print("|---|---|---|---|---|")
    for (v, m), j in sorted(E.items(), key=lambda t: (t[0][1], t[1]["final_mse"])):
        if "eps" not in j or not v.endswith(("atlas", "atlas_reg211", "atlas_zero211")):
            continue
        e = np.array(j["eps"]["D21"][1:])
        nf = f"results/noise_mlp{m}.json"
        if os.path.exists(nf):
            nz = np.array(json.load(open(nf))["D21"][1:])
            ec = np.sqrt(np.maximum(e ** 2 - nz ** 2, 0))
        else:
            ec = np.full_like(e, np.nan)
        print(f"| {v} | {m} | {np.sqrt(np.mean(e**2)):.3f} | {np.sqrt(np.mean(ec**2)):.3f} | {j['final_mse']:.2e} |")


if __name__ == "__main__":
    main()
