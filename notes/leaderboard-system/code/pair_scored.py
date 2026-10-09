# paired scored-regime comparison: python scored.py JOBDIR BASE TAG [TAG...]   (JOBDIR may be a comma list)
import re, glob, sys, numpy as np
dirs, base, tags = sys.argv[1].split(","), sys.argv[2], sys.argv[3:]
R = {}
for d in dirs:
    for f in glob.glob(d + "/log_*.txt"):
        m = re.search(r"net (\d+) (\w+): raw ([0-9.e+-]+)  C/B ([0-9.]+)  adjusted ([0-9.e+-]+)", open(f).read())
        if m: R.setdefault(m.group(2), {})[int(m.group(1))] = tuple(float(m.group(k)) for k in (3, 4, 5))
b = R[base]
for t in tags:
    x = R.get(t, {})
    for half, rng in (("all", range(100)), ("0-49", range(50)), ("50-99", range(50, 100))):
        ns = [n for n in rng if n in x and n in b]
        if not ns: continue
        out = []
        for lab, i in (("raw", 0), ("C/B", 1), ("adj", 2)):
            xb, xx = np.array([b[n][i] for n in ns]), np.array([x[n][i] for n in ns]); r = xx / xb - 1
            out.append(f"{lab} {xx.mean():.4e} vs {xb.mean():.4e} ({100*(xx.mean()/xb.mean()-1):+.2f}%, per-net {100*r.mean():+.2f} +- {100*r.std(ddof=1)/np.sqrt(len(r)):.2f}, better {int((r<0).sum())}/{len(r)})")
        print(f"{t} vs {base} [{half}, n={len(ns)}]\n   " + "\n   ".join(out))
