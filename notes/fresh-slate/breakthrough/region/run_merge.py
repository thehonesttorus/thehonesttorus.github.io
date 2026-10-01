import json, sys, time
import numpy as np
sys.path.insert(0, '../../bench')
import bench, fc
S = bench.load_set('w1024_d16')
idx = [int(a) for a in sys.argv[1].split(',')]; R = int(sys.argv[2]); sw = int(sys.argv[3])
rows = []
for i in idx:
    W = bench.weights(S, i); T = S['means'][i]
    ms = []; t0 = time.time()
    p = fc.run(W, slices=2, k4mf=True, merge=R, merge_young=2, sweeps=sw, mstats=ms)
    r = dict(mlp=i, R=R, sweeps=sw, raw=float(((p[-1] - T[-1]) ** 2).mean() - S['noise'][i]), wall=time.time() - t0,
             merge_units=float(sum(m['units'] for m in ms)), merge_wall=float(sum(m['wall'] for m in ms)), merges=len(ms),
             N_last=ms[-1]['N'] if ms else 0)
    rows.append(r); print(json.dumps(r), flush=True)
json.dump(rows, open(f'results/merge_R{R}_s{sw}_{sys.argv[1].replace(",", "")}.json', 'w'), indent=1)
