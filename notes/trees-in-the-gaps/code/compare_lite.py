import numpy as np, sys, time
from closure import closure
from lite import lite
n, L, s = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
Ws = list(np.load(f"W_n{n}_L{L}_s{s}.npy")); tr = np.load(f"truth_n{n}_L{L}_s{s}.npz")
m, v, T = tr["m"], tr["v"], int(tr["T"])
runs = {"closure": closure(Ws)}
for name, kw in {"D3+D4": dict(use=("D3","D4")), "+PP4": dict(use=("D3","D4","PP4")), "+K22": dict(use=("D3","D4","PP4","K22")),
                 "lite(all)": dict()}.items():
    t0 = time.time(); runs[name] = lite(Ws, **kw); print(name, "%.1fs" % (time.time()-t0), flush=True)
print("layer  MC@B     " + "  ".join(f"{k:>10s}" for k in runs) + "   noise")
for l in range(L):
    print(f"{l+1:3d}  {v[l].mean()/65536:.2e}  " + "  ".join(f"{np.mean((o[l]['m']-m[l])**2):10.2e}" for o in runs.values()) + f"  {v[l].mean()/T:.1e}")
for k, o in runs.items():
    e = o[L-1]["m"]-m[L-1]; print(f"{k:10s} final: bias {e.mean():+.2e} rms {np.sqrt(np.mean(e**2)):.2e}")
