import numpy as np, sys
from closure import closure
from corrected2 import corrected_closure2
n, L, s = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
Ws = list(np.load(f"W_n{n}_L{L}_s{s}.npy")); tr = np.load(f"truth_n{n}_L{L}_s{s}.npz")
m, v, T = tr["m"], tr["v"], int(tr["T"])
runs = {"closure": closure(Ws)}
for rm in [0, 1, 2, 99]:
    runs[f"r<={rm}"] = corrected_closure2(Ws, rmax=rm)
runs["r<=99 novar"] = corrected_closure2(Ws, rmax=99, var_corr=False)
print("layer  MC@B     " + "  ".join(f"{k:>11s}" for k in runs) + "   noise")
for l in range(L):
    print(f"{l+1:3d}  {v[l].mean()/65536:.2e}  " + "  ".join(f"{np.mean((o[l]['m']-m[l])**2):11.2e}" for o in runs.values()) + f"  {v[l].mean()/T:.1e}")
