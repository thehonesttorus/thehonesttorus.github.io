# Closure + a small Monte Carlo used only to fit a few coherent residual modes (the coherent part of the ideal).
import numpy as np, sys
from closure import closure, phi
from scipy.special import ndtr
n, L, s, Tc = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
Ws = list(np.load(f"W_n{n}_L{L}_s{s}.npy")); tr = np.load(f"truth_n{n}_L{L}_s{s}.npz")
m, v = tr["m"][-1], tr["v"][-1]
o = closure(Ws)[-1]; a = o["mu"]/o["sig"]; sg = o["sig"]
X = np.array([np.ones(n), sg*phi(a), sg*a*phi(a), sg*(a*a-1)*phi(a), sg*ndtr(a), o["m"]]).T
rng = np.random.default_rng(3)
base = np.mean((o["m"]-m)**2); mcb = v.mean()/65536
res = {k: [] for k in [1, 2, 4, 6]}; ora = {}
for rep in range(200):
    yhat = m + rng.standard_normal(n)*np.sqrt(v/Tc)    # MC means with Tc samples (independent across neurons approx.)
    r = o["m"] - yhat
    for k in res:
        b, *_ = np.linalg.lstsq(X[:, :k], r, rcond=None)
        res[k].append(np.mean((o["m"] - X[:, :k] @ b - m)**2))
for k in res:
    b, *_ = np.linalg.lstsq(X[:, :k], o["m"]-m, rcond=None); ora[k] = np.mean((o["m"]-X[:, :k]@b-m)**2)
print(f"n={n} L={L}: closure MSE {base:.2e} (MC@B {mcb:.2e}); calibration with {Tc} MC samples ({Tc/65374:.3f} of budget):")
for k in res: print(f"   {k} modes: MSE {np.mean(res[k]):.2e}  (oracle fit {ora[k]:.2e})")
