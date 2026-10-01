"""Team G, test T3: width universality of the age-resolution law k(a) = c n / a.
For each network: FC (exact ages), then the oracle age multiresolution at several c; reports raw MSE (truth noise
subtracted) and the noise-free truncation loss delta = mean((p_c - p_FC)^2) on the final layer.
c = 0 means: drop every age > 2 (window 2), the normaliser of the loss.
usage: python run_g.py SET MLPS CS TAG [diag]"""
import json, sys, time, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '../../fresh-slate/bench'))
import bench, fcg
S = bench.load_set(sys.argv[1])
idx = [int(a) for a in sys.argv[2].split(',')]
cs = [float(c) for c in sys.argv[3].split(',')] if sys.argv[3] != '-' else []
tag = sys.argv[4]; want_diag = len(sys.argv) > 5
out = []
for i in idx:
    W = bench.weights(S, i); T = S['means'][i]; nz = float(S['noise'][i])
    dg = [] if want_diag else None
    t0 = time.time(); p0 = fcg.run(W, slices=2, k4mf=True, agediag=dg); dt = time.time() - t0
    r = dict(set=sys.argv[1], mlp=i, c=None, raw=float(((p0[-1] - T[-1]) ** 2).mean() - nz), delta=0.0, wall=dt)
    if dg is not None: r['diag'] = dg
    out.append(r); print(json.dumps({k: v for k, v in r.items() if k != 'diag'}), flush=True)
    for c in cs:
        t0 = time.time(); p = fcg.run(W, slices=2, k4mf=True, **(dict(window=2) if c == 0 else dict(agek=c))); dt = time.time() - t0
        r = dict(set=sys.argv[1], mlp=i, c=c, raw=float(((p[-1] - T[-1]) ** 2).mean() - nz),
                 delta=float(((p[-1] - p0[-1]) ** 2).mean()), wall=dt)
        out.append(r); print(json.dumps(r), flush=True)
    json.dump(out, open(os.path.join(HERE, f'results/{tag}.json'), 'w'), indent=1, default=float)
