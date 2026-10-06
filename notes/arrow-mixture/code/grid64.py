import numpy as np, sys, time
sys.path.insert(0, "../num12")
from arrow import arrow
from closure import closure
from gac import gac
mse = lambda a, b: np.mean((a-b)**2)
seeds = [100 + i for i in range(int(sys.argv[1]))]
cfgs = [("closure",), ("gac",)] + [tuple(c.split(":")) for c in sys.argv[2:]]
res = {c: [] for c in cfgs}
for s in seeds:
    Ws = list(np.load(f"../num13/W_n64_L8_s{s}.npy")); mt = np.load(f"../num13/truth_n64_L8_s{s}.npz")["m"]
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
    print(f"seed {s}: " + "  ".join(f"{':'.join(c)} {res[c][-1]:.2e}" for c in cfgs), flush=True)
print("MEAN over seeds:")
for c in cfgs:
    print(f"  {':'.join(c):18s} {np.mean(res[c]):.3e}   ratio to closure {np.mean(res[c])/np.mean(res[cfgs[0]]):.3f}")
