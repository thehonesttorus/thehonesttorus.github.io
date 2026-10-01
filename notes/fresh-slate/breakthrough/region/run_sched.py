import json, sys, time
import numpy as np
sys.path.insert(0, '../../bench')
import bench, fc
S = bench.load_set('w1024_d16')
idx = [int(a) for a in sys.argv[1].split(',')]; R = int(sys.argv[2]); cuts = [int(c) for c in sys.argv[3].split(',')]; lo = int(sys.argv[4])
rows = []
for i in idx:
    W = bench.weights(S, i); T = S['means'][i]
    bs = []; t0 = time.time()
    p = fc.run(W, slices=2, k4mf=True, binsched=dict(cuts=cuts, lo=lo, R=R, kind='tucker'), binstats=bs)
    r = dict(mlp=i, R=R, cuts=cuts, lo=lo, raw=float(((p[-1] - T[-1]) ** 2).mean() - S['noise'][i]), wall=time.time() - t0,
             bins=[dict(N=b['N'], cut_units=b['cut_units']) for b in bs])
    rows.append(r); print(json.dumps(r), flush=True)
json.dump(rows, open(f'results/sched_R{R}_lo{lo}_c{"-".join(map(str, cuts))}_{sys.argv[1].replace(",", "")}.json', 'w'), indent=1)
