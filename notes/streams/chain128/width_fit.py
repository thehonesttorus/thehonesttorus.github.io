"""Summarise results/ladder_w*_mlp0.txt: rms over layers 1-6 (common to all widths) of the noise-corrected D21(l+1)
error per rule, and a power-law fit eps ~ n^-p extrapolated to n = 1024."""
import glob, re, sys
import numpy as np

prefix = sys.argv[1] if len(sys.argv) > 1 else "results/ladder_w"

rows = {}
for f in sorted(glob.glob(prefix + "*_mlp0.txt")):
    n = int(re.search(r"_w(\d+)_", f).group(1))
    lines = [l.split() for l in open(f) if l[:1].isdigit()]
    hdr = [l for l in open(f) if l.startswith("l ")][0].split()[2:]
    vals = np.array([[float(x) for x in l[2:]] for l in lines])
    rows[n] = (hdr, vals)
ns = sorted(rows)
hdr = rows[ns[0]][0]
print("| rule | " + " | ".join(f"n={n}" for n in ns) + " | fitted p (eps ~ n^-p) | extrapolated n=1024 |")
print("|---|" + "---|" * len(ns) + "---|---|")
for j, r in enumerate(hdr):
    if r == "true":
        continue
    e = [float(np.sqrt(np.mean(rows[n][1][1:7, j] ** 2))) for n in ns]
    p, c = np.polyfit(np.log(ns), np.log(e), 1)
    print(f"| {r} | " + " | ".join(f"{x:.3f}" for x in e) + f" | {-p:.2f} | {np.exp(c) * 1024 ** p:.3f} |")
