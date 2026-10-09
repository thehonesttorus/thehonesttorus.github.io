# paired raw comparison from mjob logs: python pair.py LOG BASE TAG [TAG...]
import re, sys, numpy as np
log, base, tags = sys.argv[1], sys.argv[2], sys.argv[3:]
R = {}
for l in open(log):
    m = re.match(r"\[(\w+)_(\d+) exit 0\] net \d+ \w+: raw ([0-9.e+-]+)", l)
    if m: R.setdefault(m.group(1), {})[int(m.group(2))] = float(m.group(3))
b = R.get(base, {})
for t in tags:
    x = R.get(t, {}); ns = sorted(set(x) & set(b))
    if not ns: print(f"{t}: no paired nets"); continue
    r = np.array([x[n] / b[n] - 1 for n in ns])
    print(f"{t:6s} vs {base}: n={len(ns):2d}  mean raw {np.mean([x[n] for n in ns]):.4e} vs {np.mean([b[n] for n in ns]):.4e}  "
          f"per-net {100*r.mean():+.2f}% +- {100*r.std(ddof=1)/np.sqrt(len(r)):.2f}  better on {(r<0).sum()}/{len(r)}")
