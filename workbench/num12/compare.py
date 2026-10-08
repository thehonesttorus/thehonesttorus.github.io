import numpy as np, sys
from closure import closure
n, L, s = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
Ws = np.load(f"W_n{n}_L{L}_s{s}.npy"); tr = np.load(f"truth_n{n}_L{L}_s{s}.npz")
m, v, T = tr["m"], tr["v"], int(tr["T"])
out = closure(list(Ws))
Tb = 65536
for l in range(L):
    e = out[l]["m"] - m[l]
    mse = (e**2).mean(); mc = v[l].mean()/Tb; noise = v[l].mean()/T
    print(f"layer {l+1:2d}: closure MSE {mse:.3e}  (truth noise {noise:.1e})  MC@B MSE {mc:.3e}  ratio closure/MC {mse/mc:6.2f}  mean err {e.mean():+.2e} rel {np.sqrt(mse)/m[l].mean():.2e}")
