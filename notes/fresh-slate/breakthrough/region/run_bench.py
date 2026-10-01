"""FC (self-contained: first-order chaos sources, exact slices, mean-field kappa4) on the bench set w1024_d16, paired with
the exact Gaussian closure.  usage: python run_bench.py MLPS "dict(...)" TAG"""
import json, sys, time
import numpy as np
sys.path.insert(0, '../../bench')
import bench, fc, gclose

S = bench.load_set('w1024_d16')
idx = [int(a) for a in sys.argv[1].split(',')] if len(sys.argv) > 1 else range(len(S['seeds']))
kw = eval(sys.argv[2]) if len(sys.argv) > 2 else dict(slices=2, k4mf=True)
rows = []
for i in idx:
    W = bench.weights(S, i); T = S['means'][i]
    t0 = time.time(); p = fc.run(W, **kw); dt = time.time() - t0
    g = gclose.predict_exact(W)
    r = dict(mlp=i, fc_raw=float(((p[-1] - T[-1]) ** 2).mean() - S['noise'][i]),
             gauss_raw=float(((g[-1] - T[-1]) ** 2).mean() - S['noise'][i]), noise=float(S['noise'][i]),
             fc_layers=((p - T) ** 2).mean(1).tolist(), wall=dt)
    rows.append(r)
    print(json.dumps({k: v for k, v in r.items() if k != 'fc_layers'}), flush=True)
fr = np.array([r['fc_raw'] for r in rows]); gr = np.array([r['gauss_raw'] for r in rows])
summ = dict(kw=str(kw), fc_raw=float(fr.mean()), fc_se=float(fr.std(ddof=1) / np.sqrt(len(fr))) if len(fr) > 1 else None,
            gauss_raw=float(gr.mean()), ratio_geo=float(np.exp(np.mean(np.log(gr / np.maximum(fr, 1e-12))))), rows=rows)
print('SUMMARY', json.dumps({k: v for k, v in summ.items() if k != 'rows'}))
tag = sys.argv[3] if len(sys.argv) > 3 else 'fc_mf'
json.dump(summ, open(f'results/bench_{tag}.json', 'w'), indent=1)
