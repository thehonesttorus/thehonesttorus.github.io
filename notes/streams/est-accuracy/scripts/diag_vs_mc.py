"""Compare the chain's per-layer pre-activation internals (lab_run EA_DUMP) with Monte-Carlo truth (mc_moments):
relative error of mean, variance, kappa3 diagonal (D3) and kappa4 diagonal (g4 = 2 dG), with the MC noise level
(from the A/B half-samples).

  python diag_vs_mc.py --dump DUMPDIR/TAG_mlpI.npz --mc results/mc/mc_I.npz
"""
import argparse
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--dump", required=True)
ap.add_argument("--mc", required=True)
a = ap.parse_args()
d = np.load(a.dump)
m = np.load(a.mc)


def rel(x, y):
    return float(np.sqrt(np.mean((x - y) ** 2) / np.mean(y ** 2)))


print("layer | mean err | var err (mean ratio) | D3 err (ratio, noise) | k4 err (ratio of means, noise)")
for l in range(16):
    out = [f"{l:2d}"]
    out.append(f"{rel(d[f'mu_{l}'], m['mean_all'][l]):.2e}")
    v = d[f"var_{l}"]; vt = m["var_all"][l]
    out.append(f"{rel(v, vt):.2e} ({v.mean() / vt.mean():.4f})")
    if f"D3_{l}" in d:
        k3 = d[f"D3_{l}"]; k3t = m["k3_all"][l]
        nz = rel(m["k3_A"][l], m["k3_B"][l]) / np.sqrt(2) / np.sqrt(2)
        out.append(f"{rel(k3, k3t):.2f} ({np.dot(k3, k3t) / np.dot(k3t, k3t):.3f}, {nz:.2f})")
    else:
        out.append("-")
    if f"g4_{l}" in d:
        g4 = d[f"g4_{l}"]; k4t = m["k4_all"][l]
        nz = rel(m["k4_A"][l], m["k4_B"][l]) / 2
        out.append(f"{rel(g4, k4t):.2f} ({g4.mean() / k4t.mean():.3f}, {nz:.2f})")
    else:
        out.append("-")
    print(" | ".join(out))
