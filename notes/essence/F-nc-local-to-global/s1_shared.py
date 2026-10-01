"""Team F test S1: old sources (age >= amin) read through ONE shared target-space basis of size k per dyadic age
bin (oracle SVD), against each source's own basis of the same size, inside costate's exact first-order co-state."""
import sys, os, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "../../fresh-slate/breakthrough/costate")); sys.path.insert(0, os.path.join(HERE, "../../fresh-slate/bench"))
import bench, costate_f as cf
name, i = sys.argv[1], int(sys.argv[2]); specs = sys.argv[3].split(",")
S = bench.load_set(name); W = bench.weights(S, i).astype(float); truth = S["means"][i]; noise = S["noise"][i]; n = W.shape[1]
for sp in specs:
    kind, amin, k = sp.split(":"); amin = int(amin)
    kk = (lambda a: None if a < amin else (min(n, int(2 * n / a)) if k == "2n/a" else int(k)))
    kw = dict(shared=kk) if kind == "shared" else (dict(rank=kk) if kind == "own" else {})
    pred = cf.predict(W, **kw)
    raw = float(((pred[-1] - truth[-1]) ** 2).mean() - noise)
    print(json.dumps(dict(mlp=i, spec=sp, raw=raw)), flush=True)
