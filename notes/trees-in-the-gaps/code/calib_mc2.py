# Residue calibration with the radial trick: E h(x) = E|x| * E h(x/|x|) exactly (h is 1-homogeneous).
import numpy as np, sys
from math import lgamma, exp, sqrt
from closure import closure
n, L, s = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
Ws = list(np.load(f"W_n{n}_L{L}_s{s}.npy")); tr = np.load(f"truth_n{n}_L{L}_s{s}.npz"); m = tr["m"][-1]; v = tr["v"][-1]
mc = closure(Ws)[-1]["m"]; Wf = [W.T.astype(np.float32) for W in Ws]
Enorm = sqrt(2)*exp(lgamma((n+1)/2)-lgamma(n/2))
rng = np.random.default_rng(11)
print(f"n={n} L={L}: closure MSE {np.mean((mc-m)**2):.2e}, oracle-scale {np.mean((mc*(m@mc)/(mc@mc)-m)**2):.2e}, MC@B {v.mean()/65536:.2e}")
for Tc in [250, 1000, 4000]:
    r1, r2 = [], []
    for rep in range(40):
        x = rng.standard_normal((Tc, n)).astype(np.float32)
        xr = x/np.linalg.norm(x, axis=1, keepdims=True)*Enorm
        for X, out in [(x, r1), (xr, r2)]:
            h = X
            for W in Wf: h = np.maximum(h @ W, 0)
            y = h.mean(0); c = mc @ (mc - y)/(mc @ mc)
            out.append(np.mean(((1-c)*mc - m)**2))
    print(f"   Tc={Tc:5d} ({Tc/65374*100:.1f}%): plain {np.mean(r1):.2e}   radial {np.mean(r2):.2e}")
# per-neuron MC variance reduction from the radial trick
x = rng.standard_normal((20000, n)).astype(np.float32); xr = x/np.linalg.norm(x, axis=1, keepdims=True)*Enorm
hs = []
for X in [x, xr]:
    h = X
    for W in Wf: h = np.maximum(h @ W, 0)
    hs.append(h)
print("   per-neuron variance ratio radial/plain: %.3f; variance of projection on m: ratio %.3f" % (hs[1].var(0).mean()/hs[0].var(0).mean(), (hs[1]@mc).var()/(hs[0]@mc).var()))
