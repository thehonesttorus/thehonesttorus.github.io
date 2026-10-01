import json, sys, time
import numpy as np
sys.path.insert(0, '../../bench')
import bench, fc
S = bench.load_set('w1024_d16')
idx = [int(a) for a in sys.argv[1].split(',')]; R = int(sys.argv[2]); kind = sys.argv[3]; m = int(sys.argv[4]) if len(sys.argv) > 4 else 9
rows = []
for i in idx:
    W = bench.weights(S, i); T = S['means'][i]
    bs = []; t0 = time.time()
    p = fc.run(W, slices=2, k4mf=True, bincut=dict(m=m, lo=5, hi=8, R=R, kind=kind), binstats=bs)
    r = dict(mlp=i, R=R, kind=kind, cut=m, raw=float(((p[-1] - T[-1]) ** 2).mean() - S['noise'][i]), wall=time.time() - t0,
             cut_stats=bs[0], per_layer=bs[1:])
    rows.append(r)
    print(json.dumps(dict(mlp=i, R=R, kind=kind, raw=r['raw'], cut_units=bs[0]['cut_units'], cut_wall=bs[0]['wall'],
                          eps_bin=[round(x['eps_bin'], 3) for x in bs[1:]], bin_share=[round(x['bin_share'], 3) for x in bs[1:]],
                          eps_total=[round(x['eps_total'], 3) for x in bs[1:]])), flush=True)
json.dump(rows, open(f'results/bin_{kind}_R{R}_m{m}_{sys.argv[1].replace(",", "")}.json', 'w'), indent=1)
