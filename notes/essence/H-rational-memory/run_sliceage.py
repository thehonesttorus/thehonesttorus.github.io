"""Score-level test: FC with slice legs only for ages <= K, on the w1024_d16 bench (paired with full FC from region).
usage: python3 run_sliceage.py MLPS K TAG"""
import json, sys, time
import numpy as np
sys.path.insert(0, '../../fresh-slate/bench')
import bench, fc_age
S = bench.load_set('w1024_d16')
mlps = [int(a) for a in sys.argv[1].split(',')]
K = None if sys.argv[2] == 'all' else int(sys.argv[2])
rows = []
for i in mlps:
    W = bench.weights(S, i); T = S['means'][i]
    cost = []
    t0 = time.time(); p = fc_age.run(W, slices=2, k4mf=True, slice_age=K, cost=cost); dt = time.time() - t0
    c = np.array([x[2:] for x in cost])
    r = dict(mlp=i, K=K, raw=float(((p[-1] - T[-1]) ** 2).mean() - S['noise'][i]), wick_products=float(c[:, 0].sum()),
             slice_products=float(c[:, 1].sum()), wall=round(dt, 1))
    rows.append(r); print(json.dumps(r), flush=True)
json.dump(rows, open(f'results/sliceage_{sys.argv[3]}.json', 'w'), indent=1)
