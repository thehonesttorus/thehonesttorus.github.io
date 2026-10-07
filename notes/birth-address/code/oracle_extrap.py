# Oracle attribution: noise-free extrapolation MSE(N) = a + b/N from the half-sample and full-sample oracle runs.
#   python oracle_extrap.py oracle_off0.txt [oracle_off1.txt ...]
# Lines "net K ORACLE:h0|h1|full: raw X" (several files of one network may be given; later lines win).
import re, sys
r = {}; net = None
for p in sys.argv[1:]:
    for l in open(p):
        m = re.match(r"net (\d+) (\S+): raw (\S+)", l)
        if m:
            r[m[2]] = float(m[3]); net = m[1]
base = r.get("chain", r.get("none:none"))
print(f"net {net}: chain {base:.4e}")
orcs = sorted({k.rsplit(":", 1)[0] for k in r if k.endswith(":full") and not k.startswith("none")}, key=lambda o: (o.count("+"), o))
for o in orcs:
    h = [r[f"{o}:h{k}"] for k in (0, 1) if f"{o}:h{k}" in r]
    full = r[f"{o}:full"]
    a = 2 * full - sum(h) / len(h) if h else float("nan")
    print(f"  oracle {o:28s}: halves {' '.join(f'{x:.4e}' for x in h):23s} | full {full:.4e} | noise-free a = {a:.4e} "
          f"({100 * (a / base - 1):+6.1f}% vs chain; noise at full N {100 * (full - a) / base:+5.1f}%)")
