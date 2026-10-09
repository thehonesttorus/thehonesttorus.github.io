# Directional heat defect of the adopted chain (note XLII, the F1 gating experiment):
#   V60_F64=1 V33_SAT=nan python heatdef.py NET I0 I1 [H=0.5] [SEED=0]
# For unit directions v_r (r = I0..I1-1 of a list drawn from SEED + NET), with E(m, Sigma) the chain on N(m, Sigma):
#   d_r = [E(0, I + h v v^T) - E(0, I - h v v^T)] / (2h) - [E(h v, I) + E(-h v, I) - 2 E(0, I)] / (2 h^2)
# (the heat defect tr(D P_v), P_v = v v^T, its covariance part acting as the control variate of the mean part).
# Writes $OUT/heatdef_{NET}_{I0}.npz: d (R x 16 x n), cov part, mean part, E0, the directions' indices, h.
import os, sys, time, importlib.util, numpy as np, flopscope as flops
from whestbench import MLP
net, i0, i1 = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
h = float(sys.argv[4]) if len(sys.argv) > 4 else 0.5
seed = int(sys.argv[5]) if len(sys.argv) > 5 else 0
n = 1024
V = np.random.default_rng(7919 * (seed + 1) + net).standard_normal((i1, n))
V /= np.linalg.norm(V, axis=1, keepdims=True)
spec = importlib.util.spec_from_file_location("estf", "estimator_final_v56.py")
mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
Wcol = np.load(f"../official/W_off{net}.npy")
mlp = MLP(width=n, depth=16, weights=[np.ascontiguousarray(W.T).astype(np.float32) for W in Wcol])
est = mod.Estimator()


def run(m=None, cov=None):
    est.in_mean, est.in_cov = m, cov
    with flops.BudgetContext(flop_budget=2**46, wall_time_limit_s=6000.0, quiet=True):
        return np.asarray(est.predict(mlp, 2**41), dtype=np.float64)


t0 = time.time()
E0 = run()
d, cp, mp = [], [], []
for r in range(i0, i1):
    v = V[r]
    c = (run(cov=(h, v)) - run(cov=(-h, v))) / (2 * h)
    m = (run(m=h * v) + run(m=-h * v) - 2 * E0) / (2 * h * h)
    d.append(c - m); cp.append(c); mp.append(m)
    print(f"net {net}: direction {r} done at {time.time() - t0:.0f}s", flush=True)
np.savez(f"{os.environ.get('OUT', '.')}/heatdef_{net}_{i0}.npz", d=np.array(d), cov=np.array(cp), mean=np.array(mp),
         E0=E0, idx=np.arange(i0, i1), h=h)
print(f"net {net}: directions {i0}-{i1 - 1} in {time.time() - t0:.0f}s", flush=True)
