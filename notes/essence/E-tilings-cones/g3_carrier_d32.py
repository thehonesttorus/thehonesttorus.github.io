"""G3 test: team D's causal age-multiresolution carrier (fcs.run(causal=c), read-only import) on the smoke shape
w256_d32, against full FC. Raw at the final layer (32) and at layer 16 of the same run (truth: bench all-layer means).
usage: python g3_carrier_d32.py MLPS "dict(...)" TAG"""
import json, sys, time, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '../../fresh-slate/bench')); sys.path.insert(0, os.path.join(HERE, '../D-elimination-sparsification'))
import bench, fcs
S = bench.load_set('w256_d32')
idx = [int(a) for a in sys.argv[1].split(',')]; kw = eval(sys.argv[2]); tag = sys.argv[3]
out = []
for i in idx:
    W = bench.weights(S, i); T = S['means'][i]
    t0 = time.time(); p = fcs.run(W, slices=2, k4mf=True, **kw); dt = time.time() - t0
    per = ((p - T) ** 2).mean(1)
    noise_l = S['noise'][i] * T.var(1) / T[-1].var() if False else None
    r = dict(mlp=i, raw32=float(per[-1] - S['noise'][i]), mse16=float(per[15]), per_layer=per.tolist(), wall=dt, kw=str(kw))
    out.append(r); print(json.dumps({k: v for k, v in r.items() if k != 'per_layer'}), flush=True)
os.makedirs(os.path.join(HERE, 'results'), exist_ok=True)
json.dump(out, open(os.path.join(HERE, f'results/g3_{tag}.json'), 'w'), indent=1)
