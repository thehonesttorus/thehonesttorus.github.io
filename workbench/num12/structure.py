import numpy as np, sys
from closure import closure, phi
from scipy.special import ndtr
n, L, s = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
Ws = np.load(f"W_n{n}_L{L}_s{s}.npy"); tr = np.load(f"truth_n{n}_L{L}_s{s}.npz")
m = tr["m"]; out = closure(list(Ws))
for l in [1, L//2, L-1]:
    o = out[l]; e = o["m"] - m[l]; a = o["mu"]/o["sig"]; sg = o["sig"]
    feats = {"1": np.ones(n), "sig*phi(a)": sg*phi(a), "sig*a*phi": sg*a*phi(a), "sig*(a2-1)phi": sg*(a*a-1)*phi(a), "m": o["m"], "sig*Phi": sg*ndtr(a)}
    X = np.array(list(feats.values())).T
    for k in [1, 2, 4, 6]:
        Xk = X[:, :k]; c, *_ = np.linalg.lstsq(Xk, e, rcond=None)
        r = e - Xk @ c
        print(f"layer {l+1}: features {list(feats)[:k]}: explained {1-(r**2).mean()/(e**2).mean():.3f}")
    print(f"   a range [{a.min():.2f},{a.max():.2f}] mean {a.mean():.2f}")
