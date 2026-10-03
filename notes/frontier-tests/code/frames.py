# Structured angular quadrature at width 1024: signed-Hadamard frames + antipodes + exact radius, vs i.i.d. Gaussian,
# at matched forward passes (~10% of the budget = 6144 passes). Several independent seeds -> MSE vs 1e9 truth.
import numpy as np, sys
from math import lgamma, exp, sqrt
from scipy.linalg import hadamard
i = int(sys.argv[1]) if len(sys.argv) > 1 else 0
Wcol = np.load(f"../official/W_off{i}.npy"); mt = np.load(f"../official/truth_off{i}.npz")["m"][-1].astype(np.float64)
Wf = [np.ascontiguousarray(W.T).astype(np.float32) for W in Wcol]; n = 1024
Rbar = sqrt(2)*exp(lgamma((n+1)/2) - lgamma(n/2)); H = (hadamard(n)/np.sqrt(n)).astype(np.float32)
def run(X):
    h = X
    for W in Wf: h = np.maximum(h @ W, 0)
    return h.astype(np.float64).mean(0)
res = {"iid Gaussian": [], "iid Gaussian, antithetic": [], "signed-Hadamard frames + antipodes, exact radius": []}
for seed in range(6):
    rng = np.random.default_rng(100 + seed)
    X = rng.standard_normal((6144, n)).astype(np.float32); res["iid Gaussian"].append(np.mean((run(X) - mt)**2))
    X = rng.standard_normal((3072, n)).astype(np.float32); res["iid Gaussian, antithetic"].append(np.mean((run(np.concatenate([X, -X])) - mt)**2))
    F = np.concatenate([H * rng.choice([-1, 1], n).astype(np.float32)[None, :] @ np.linalg.qr(rng.standard_normal((n, n)))[0].astype(np.float32) for _ in range(3)])
    res["signed-Hadamard frames + antipodes, exact radius"].append(np.mean((run(np.concatenate([F, -F])*Rbar) - mt)**2))
for k, v in res.items(): print(f"{k:50s} MSE {np.mean(v):.3e} (+-{np.std(v)/np.sqrt(len(v)):.1e})")
