"""Separator test on the bench: final-layer raw MSE (4 MLPs of w128_d16) for history separators."""
import sys, json, numpy as np
sys.path.insert(0, '.'); sys.path.insert(0, '/home/user/thehonesttorus.github.io/notes/fresh-slate/bench')
import bench
from mkv import mkv
S = bench.load_set(sys.argv[1]); idx = range(int(sys.argv[2]))
V = {'diag w1': dict(w=1), 'rank8 w1': dict(w=1, rank=8), 'rank32 w1': dict(w=1, rank=32),
     'diag w2': dict(w=2), 'rank16 w2': dict(w=2, rank=16), 'full': dict(w=16)}
res = {}
for i in idx:
    W = bench.weights(S, i); t = S['means'][i]; nz = S['noise'][i] if np.ndim(S['noise']) else S['noise']
    for k, kw in V.items():
        e = ((mkv(W, var21='full', **kw) - t) ** 2).mean(1)
        res.setdefault(k, []).append(float(e[-1] - nz))
        print(i, k, f"{e[-1]-nz:.2e}", flush=True)
for k, r in res.items(): print(f"== {k:10s} mean {np.mean(r):.2e}")
json.dump(res, open(f'results/sep_{sys.argv[1]}.json', 'w'))
