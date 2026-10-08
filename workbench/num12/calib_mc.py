# One-scalar ("residue") calibration of the closure with *actual* correlated Monte Carlo samples.
import numpy as np, sys
from closure import closure
n, L, s = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
Ws = list(np.load(f"W_n{n}_L{L}_s{s}.npy")); tr = np.load(f"truth_n{n}_L{L}_s{s}.npz"); m = tr["m"][-1]; v = tr["v"][-1]
mc = closure(Ws)[-1]["m"]; Wf = [W.T.astype(np.float32) for W in Ws]
rng = np.random.default_rng(11)
print(f"n={n} L={L}: closure MSE {np.mean((mc-m)**2):.2e}, oracle-scale {np.mean((mc*(m@mc)/(mc@mc)-m)**2):.2e}, MC@B {v.mean()/65536:.2e}")
for Tc in [250, 1000, 4000]:
    res = []
    for rep in range(40):
        h = rng.standard_normal((Tc, n)).astype(np.float32)
        for W in Wf: h = np.maximum(h @ W, 0)
        y = h.mean(0)
        c = mc @ (mc - y)/(mc @ mc)
        res.append(np.mean(((1-c)*mc - m)**2))
    print(f"   Tc={Tc:5d} ({Tc/65374*100:.1f}% budget): calibrated MSE {np.mean(res):.2e} +- {np.std(res)/np.sqrt(len(res)):.1e}")
