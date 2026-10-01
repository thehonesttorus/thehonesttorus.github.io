import sys, numpy as np, json
from common import bench
import tc
sets = sys.argv[1].split(","); variants = dict(
    closureK2=dict(inject=False), tc_mean_only=dict(cov_inject=False), tc=dict(), tc_chi=dict(chi=True))
res = {}
for name in sets:
    S = bench.load_set(name)
    for i in range(len(S["seeds"])):
        W = bench.weights(S, i); T = S["means"][i]; nz = S["noise"][i]; row = []
        for k, kw in variants.items():
            p = tc.predict(W, **kw); r = ((p[-1] - T[-1]) ** 2).mean() - nz
            res.setdefault((name, k), []).append(r); row.append(f"{k} {r:.3e}")
        print(name, i, " | ".join(row), flush=True)
for (name, k), v in res.items(): print(name, k, f"{np.mean(v):.3e}")
