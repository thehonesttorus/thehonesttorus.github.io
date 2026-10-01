"""Team G: score-level comparison of age-resolution profiles against the uniform law k = c n / a at the SAME old-atom
count. The uniform loss at the matching c is interpolated log-linearly in c from the T3 sweep (the loss is exponential in c).
usage: python compare_prof.py PROF_TAG N T3_TAG [T3_TAG ...]"""
import json, sys
import numpy as np


def make(pr):
    t = pr.split(':')
    if t[0] == 'lin':
        c, g, fm = map(float, t[1:]); return lambda a, n, L: n * max(c / a - g, fm)
    if t[0] == 'pow':
        c, p = map(float, t[1:]); return lambda a, n, L: n * min(1.0, c / a ** p)
    raise ValueError(pr)


def atoms(fn, n, L):  # same count as run_prof.atoms
    s = sum(min(n, max(1, round(fn(a, n, L)))) for t in range(L) for a in range(3, t + 2))
    s2 = sum(min(n, round(2 * n / a)) for t in range(L) for a in range(3, t + 2))
    return s / s2


ptag, n, ttags = sys.argv[1], int(sys.argv[2]), sys.argv[3:]
T = {}
for t in ttags:
    for r in json.load(open(f'results/{t}.json')):
        T.setdefault(r['mlp'], {})[r['c']] = r['delta']
cg = np.linspace(0.5, 3, 251)
ag = np.array([atoms(lambda a, n_, L, c=c: n_ * min(1.0, c / a), n, 16) for c in cg])
out = {}
print('==', ptag)
for r in json.load(open(f'results/{ptag}.json')):
    if r['prof'] is None:
        continue
    d = T[r['mlp']]; cs = sorted(c for c in d if c not in (None, 0.0)); ls = np.log([d[c] for c in cs])
    cu = float(np.interp(r['old_atoms_vs_c2'], ag, cg)); lu = float(np.exp(np.interp(cu, cs, ls)))
    out.setdefault(r['prof'], []).append(r['delta'] / lu)
    print(f" mlp {r['mlp']} {r['prof']:22s} atoms {r['old_atoms_vs_c2']:.3f} loss {r['delta']:.3g} | uniform c = {cu:.2f}: "
          f"loss {lu:.3g} | ratio {r['delta'] / lu:.2f}")
for p, v in out.items():
    print(f'  {p:22s} loss / uniform loss at equal atoms: mean {np.mean(v):.2f} (range {min(v):.2f}-{max(v):.2f})')
