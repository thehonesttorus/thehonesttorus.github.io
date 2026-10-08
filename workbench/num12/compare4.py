import numpy as np, sys, time
from closure import closure
from corrected3 import tree_closure
n, L, s = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
Ws = list(np.load(f"W_n{n}_L{L}_s{s}.npy")); tr = np.load(f"truth_n{n}_L{L}_s{s}.npz")
m, v, T = tr["m"], tr["v"], int(tr["T"])
runs = {"closure": closure(Ws)}
t0 = time.time()
runs["1pt r<=2"] = tree_closure(Ws, two_point=False)
runs["1pt+2pt r<=2"] = tree_closure(Ws)
runs["1pt+2pt r<=0"] = tree_closure(Ws, rmax=0)
print("time %.1fs" % (time.time()-t0))
print("layer  MC@B     " + "  ".join(f"{k:>13s}" for k in runs) + "   noise")
for l in range(L):
    print(f"{l+1:3d}  {v[l].mean()/65536:.2e}  " + "  ".join(f"{np.mean((o[l]['m']-m[l])**2):13.2e}" for o in runs.values()) + f"  {v[l].mean()/T:.1e}")
