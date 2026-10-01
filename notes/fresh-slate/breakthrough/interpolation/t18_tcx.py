import sys, numpy as np
from common import bench
import tc
for name in sys.argv[1].split(","):
    S = bench.load_set(name); R = {}
    for i in range(len(S["seeds"])):
        W = bench.weights(S, i); T = S["means"][i]; nz = S["noise"][i]
        for lab, kw in [("TC+chi", dict(chi=True)), ("TCX+chi", dict(chi=True, x_channel=True))]:
            p = tc.predict(W, **kw); r = ((p[-1] - T[-1]) ** 2).mean() - nz; R.setdefault(lab, []).append(r)
            R.setdefault(lab + " all-layer", []).append(((p - T) ** 2).mean())
        print(name, i, {k: f"{v[-1]:.3e}" for k, v in R.items()}, flush=True)
    print(name, {k: f"{np.mean(v):.3e}" for k, v in R.items()}, flush=True)
