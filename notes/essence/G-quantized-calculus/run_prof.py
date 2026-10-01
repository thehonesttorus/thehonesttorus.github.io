"""Team G: score-level test of age-resolution profiles beyond the uniform c n / a law.
Profile 'lin:c:g:fmin' reads a source of age a >= 3 at k(a) = n max(c/a - g, fmin); 'pow:c:p' at k = n c / a^p.
Reports raw (truth noise subtracted), the noise-free loss against FC, and the old-atom count relative to c = 2 uniform.
usage: python run_prof.py SET MLPS PROF1,PROF2 TAG"""
import json, sys, time, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '../../fresh-slate/bench'))
import bench, fcg

def make(pr):
    t = pr.split(':')
    if t[0] == 'lin':
        c, g, fm = map(float, t[1:]); return lambda a, n, L: n * max(c / a - g, fm)
    if t[0] == 'pow':
        c, p = map(float, t[1:]); return lambda a, n, L: n * min(1.0, c / a ** p)
    raise ValueError(pr)

def atoms(fn, n, L):  # sum over targets t and ages a >= 3 of k(a), relative to c = 2 uniform
    s = sum(min(n, max(1, round(fn(a, n, L)))) for t in range(L) for a in range(3, t + 2))
    s2 = sum(min(n, round(2 * n / a)) for t in range(L) for a in range(3, t + 2))
    return s / s2

S = bench.load_set(sys.argv[1]); idx = [int(a) for a in sys.argv[2].split(',')]
profs = sys.argv[3].split(','); tag = sys.argv[4]
out = []
for i in idx:
    W = bench.weights(S, i); T = S['means'][i]; nz = float(S['noise'][i]); L, n, _ = W.shape
    p0 = fcg.run(W, slices=2, k4mf=True)
    out.append(dict(set=sys.argv[1], mlp=i, prof=None, raw=float(((p0[-1] - T[-1]) ** 2).mean() - nz), delta=0.0, rel=1.0))
    print(json.dumps(out[-1]), flush=True)
    for pr in profs:
        fn = make(pr); t0 = time.time()
        p = fcg.run(W, slices=2, k4mf=True, agefn=fn)
        out.append(dict(set=sys.argv[1], mlp=i, prof=pr, raw=float(((p[-1] - T[-1]) ** 2).mean() - nz),
                        delta=float(((p[-1] - p0[-1]) ** 2).mean()), old_atoms_vs_c2=atoms(fn, n, L), wall=time.time() - t0))
        print(json.dumps(out[-1]), flush=True)
    json.dump(out, open(os.path.join(HERE, f'results/{tag}.json'), 'w'), indent=1)
