import numpy as np, sys, time
from closure import closure
from diagsrc import diag_tlp
n, L, s = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
Ws = list(np.load(f"W_n{n}_L{L}_s{s}.npy")); tr = np.load(f"truth_n{n}_L{L}_s{s}.npz")
m, v, T = tr["m"], tr["v"], int(tr["T"])
runs = {"closure": closure(Ws)}; mats = {"closure": 2*L-1}
for name, kw in {"d r0": dict(rmax=0), "d r1": dict(rmax=1), "d r2": dict(rmax=2), "d r4": dict(rmax=4), "d r1 no22": dict(rmax=1, k22=False)}.items():
    cnt = [0]; t0 = time.time(); runs[name] = diag_tlp(Ws, counter=cnt, **kw); mats[name] = cnt[0]
    print(name, "%.1fs" % (time.time()-t0), "n^3 products:", cnt[0], flush=True)
print("layer  MC@B     " + "  ".join(f"{k:>11s}" for k in runs) + "   noise")
for l in range(L):
    print(f"{l+1:3d}  {v[l].mean()/65536:.2e}  " + "  ".join(f"{np.mean((o[l]['m']-m[l])**2):11.2e}" for o in runs.values()) + f"  {v[l].mean()/T:.1e}")
B = 65374*(2*L*n*n)
for k, o in runs.items():
    e = o[L-1]["m"]-m[L-1]; frac = mats[k]*2*n**3/B
    print(f"{k:11s} final MSE {np.mean(e**2):.2e} bias {e.mean():+.2e}  n^3-products {mats[k]}  budget frac at this n {frac:.3f}")
