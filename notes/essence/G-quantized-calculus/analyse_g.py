"""Team G: T3 tables (raw, noise-free loss delta(c), relative loss delta(c)/delta(window 2)) and the per-source spectral
diagnostic: per-age energy share eta_a, participation ratio, and the resolution k_eps(a) at which the D21 contribution of a
source of age a loses a fraction eps of its energy (log-linear interpolation on the grid).
usage: python analyse_g.py tag1 tag2 ..."""
import json, sys, numpy as np
from collections import defaultdict

def keps(err, eps):
    fs = sorted(float(f) for f in err); e = [err[f] if f in err else err[str(f)] for f in fs]
    e = [max(x, 1e-30) for x in e]
    for i in range(len(fs)):
        if e[i] <= eps:
            if i == 0: return fs[0]
            # interpolate log e linearly in f
            f0, f1, l0, l1 = fs[i - 1], fs[i], np.log(e[i - 1]), np.log(e[i])
            return f0 + (np.log(eps) - l0) * (f1 - f0) / (l1 - l0)
    return 1.0

for tag in sys.argv[1:]:
    R = json.load(open(f'results/{tag}.json'))
    print(f'== {tag}')
    by = defaultdict(dict)
    for r in R:
        by[r['mlp']][r['c']] = r
    for m, d in by.items():
        base = d[None]['raw']; d0 = d.get(0.0, {}).get('delta')
        s = f"mlp {m}: FC raw {base:.3g}"
        if d0: s += f" | drop-old loss {d0:.3g}"
        for c in sorted(k for k in d if k not in (None, 0.0)):
            s += f" | c={c}: raw {d[c]['raw']:.3g}, loss {d[c]['delta']:.2g}" + (f" ({d[c]['delta'] / d0:.3f})" if d0 else '') + f", loss/FCraw {d[c]['delta'] / base:.3f}"
        print(s)
    # spectral diagnostic
    for r in R:
        if 'diag' not in r: continue
        dg = r['diag']; n = int(tag.split('w')[-1].rstrip('ab')) if 'w' in tag else None
        D2 = {x['l']: x['D2'] for x in dg if 'D2' in x}
        rows = [x for x in dg if 'age' in x]
        agg = defaultdict(lambda: defaultdict(list))
        for x in rows:
            for eps in (0.1, 0.03, 0.01):
                agg[x['age']][f'k{eps}'].append(keps(x['err'], eps))
            agg[x['age']]['pr2'].append(x['pr2'])
            agg[x['age']]['share'].append(x['Ds2'] / D2[x['l']])
        print(f" mlp {r['mlp']} per-age: a, n_src, pr2/n, a*pr2/n, energy share (mean over targets), k_eps/n and a*k_eps/n for eps=0.1,0.03,0.01")
        for a in sorted(agg):
            g = agg[a]; pr = np.mean(g['pr2'])
            # pr2 is in absolute units; normalise by n inferred from grid via share rows
            print(f"  a={a:2d} N={len(g['pr2']):2d} pr2={pr:.4g} share={np.mean(g['share']):.4f} " +
                  " ".join(f"k{e}={np.mean(g[f'k{e}']):.4f} a*k={a * np.mean(g[f'k{e}']):.3f}" for e in (0.1, 0.03, 0.01)))
