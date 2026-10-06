# Arrow-mixture grid at width n: num12 nets (n=256/1024, 1e7-1.6e7 truth) or official nets (n=1024, 1e9 truth; tag 'offK').
import numpy as np, sys, time
sys.path.insert(0, "../num12")
from arrow import arrow
from closure import closure
from gac import gac
mse = lambda a, b: np.mean((a-b)**2)
n, L = int(sys.argv[1]), int(sys.argv[2]); tags = sys.argv[3].split(",")
cfgs = [("closure",), ("gac",)] + [tuple(c.split(":")) for c in sys.argv[4:]]
res = {c: [] for c in cfgs}
for tag in tags:
    if tag.startswith("off"):
        Ws = list(np.load(f"../official/W_off{tag[3:]}.npy").astype(np.float64)); mt = np.load(f"../official/truth_off{tag[3:]}.npz")["m"]
    else:
        Ws = list(np.load(f"../num12/W_n{n}_L{L}_s{tag}.npy")); mt = np.load(f"../num12/truth_n{n}_L{L}_s{tag}.npz")["m"]
    for c in cfgs:
        t0 = time.time()
        if c[0] == "closure": out = np.array([d["m"] for d in closure(Ws)])
        elif c[0] == "gac": out = np.array([d["m"] for d in gac(Ws)])
        else:
            base, k, q, w = c[0], int(c[1]), int(c[2]), int(c[3])
            lev = int(c[4]) if len(c) > 4 and c[4] != "-" else None
            gn = len(c) > 5 and c[5] == "g"
            out = arrow(Ws, k=k, q=q, w=w, base=base, maxM=400000, lev=lev, gain=gn)
        res[c].append(mse(out[-1], mt[-1]))
        print(f"  n={n} {tag} {':'.join(c)} {res[c][-1]:.3e}  ({time.time()-t0:.0f}s)", flush=True)
print(f"MEAN over {tags} at n={n}:")
for c in cfgs:
    print(f"  {':'.join(c):22s} {np.mean(res[c]):.3e}   ratio to closure {np.mean(res[c])/np.mean(res[cfgs[0]]):.3f}")
