# Oracle attribution: noise-free extrapolation MSE(N) = a + b/N from the half-sample and full-sample oracle runs.
#   python oracle_extrap.py oracle_off0.txt [oracle_off1.txt ...]
import re, sys
for p in sys.argv[1:]:
    r = {}
    for l in open(p):
        m = re.match(r"net (\d+) (\S+): raw (\S+)", l)
        if m: r[m[2]] = float(m[3]); net = m[1]
    base = r.get("chain", r.get("none:none"))
    print(f"net {net}: chain {base:.4e}" + (f" (cold rerun {r['none:none']:.4e})" if "none:none" in r else ""))
    for o in ("D3", "D21", "D3+D21", "G4", "D3+D21+G4"):
        if f"{o}:full" not in r:
            continue
        h = [r[f"{o}:h{k}"] for k in (0, 1) if f"{o}:h{k}" in r]
        full = r[f"{o}:full"]
        a = 2 * full - sum(h) / len(h) if h else float("nan")
        print(f"  oracle {o:9s}: half-samples {' '.join(f'{x:.4e}' for x in h)} | full {full:.4e} | noise-free a = {a:.4e} "
              f"({100 * (a / base - 1):+.1f}% vs chain; noise term at full N {100 * (full - a) / base:+.1f}%)")
