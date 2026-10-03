import numpy as np, sys, time
from closure import closure
from corrected3 import tree_closure
n, L, s = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
cfgs = eval(sys.argv[4])
Ws = list(np.load(f"W_n{n}_L{L}_s{s}.npy")); tr = np.load(f"truth_n{n}_L{L}_s{s}.npz")
m, v, T = tr["m"], tr["v"], int(tr["T"])
runs = {"closure": closure(Ws)}
for name, kw in cfgs.items():
    t0 = time.time(); runs[name] = tree_closure(Ws, **kw); print(name, "%.1fs" % (time.time()-t0), flush=True)
print("layer  MC@B     " + "  ".join(f"{k:>11s}" for k in runs) + "   noise")
for l in range(L):
    print(f"{l+1:3d}  {v[l].mean()/65536:.2e}  " + "  ".join(f"{np.mean((o[l]['m']-m[l])**2):11.2e}" for o in runs.values()) + f"  {v[l].mean()/T:.1e}")
o = runs[list(cfgs)[-1]]; e = o[L-1]["m"]-m[L-1]
print("last cfg final layer: mean err %+.2e, rms %.2e, rms after removing mean %.2e" % (e.mean(), np.sqrt(np.mean(e**2)), e.std()))
