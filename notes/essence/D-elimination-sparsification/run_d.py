"""Team D tests: FC with diagnostics of the old-atom Gram structure, or with unbiased Poisson atom resampling.
usage: python run_d.py MLPS "dict(...)" TAG"""
import json, sys, time, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '../../fresh-slate/bench'))
import bench, fcs
S = bench.load_set('w1024_d16')
idx = [int(a) for a in sys.argv[1].split(',')]
kw = eval(sys.argv[2]); tag = sys.argv[3]
out = []
for i in idx:
    W = bench.weights(S, i); T = S['means'][i]
    d = [] if kw.pop('want_diag', False) or 'diag' in kw else None
    if d is not None: kw['diag'] = d
    t0 = time.time(); p = fcs.run(W, slices=2, k4mf=True, **kw); dt = time.time() - t0
    r = dict(mlp=i, raw=float(((p[-1] - T[-1]) ** 2).mean() - S['noise'][i]), wall=dt, kw=str({k: v for k, v in kw.items() if k != 'diag'}))
    if d is not None:
        r['diag'] = d
    out.append(r)
    print(json.dumps({k: v for k, v in r.items() if k != 'diag'}), flush=True)
os.makedirs(os.path.join(HERE, 'results'), exist_ok=True)
json.dump(out, open(os.path.join(HERE, f'results/{tag}.json'), 'w'), indent=1, default=float)
