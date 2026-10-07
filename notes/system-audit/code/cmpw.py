# paired comparison of run_v29w outputs: python cmpw.py BASE.txt VAR.txt [VAR2.txt ...]
import re, math, sys
def load(p):
    d = {}
    for l in open(p):
        m = re.match(r"net (\d+) \S+: raw (\S+)\s+C/B (\S+)\s+adjusted (\S+).*?residual (\S+)s\s+first-call C/B (\S+) residual (\S+)s", l)
        if m: d[int(m[1])] = tuple(float(x) for x in m.groups()[1:])
    return d
b = load(sys.argv[1])
for p in sys.argv[2:]:
    d = load(p); ks = sorted(set(b) & set(d))
    r = [d[k][0] / b[k][0] - 1 for k in ks]; m = sum(r) / len(r)
    s = math.sqrt(sum((x - m) ** 2 for x in r) / max(1, len(r) - 1) / len(r))
    a = [d[k][2] / b[k][2] - 1 for k in ks]
    print(f"{p.split('/')[-1]:>14} n={len(ks)} raw {100*m:+.3f}% +- {100*s:.3f} better {sum(x < 0 for x in r)}/{len(r)} | "
          f"C/B {sum(d[k][1] for k in ks)/len(ks):.4f} vs {sum(b[k][1] for k in ks)/len(ks):.4f} | adj {100*sum(a)/len(a):+.2f}% | "
          f"resid {sum(d[k][3] for k in ks)/len(ks):.3f} vs {sum(b[k][3] for k in ks)/len(ks):.3f} | first {sum(d[k][4] for k in ks)/len(ks):.4f} r {sum(d[k][5] for k in ks)/len(ks):.3f}")
