# Evaluate the adopted chain on shifted input laws N(m_k, s_k I) (localization / conditional-law experiments):
#   python chainat.py NET SPEC.npz [TAG]
# SPEC.npz holds M (K x n input means) and optionally s (K scales; default 1). By homogeneity the chain is run on
# N(m_k / sqrt(s_k), I) and its prediction multiplied by sqrt(s_k). Writes $OUT/chainat_{TAG}_{NET}.npz with
# out (K x 16 x n predicted means per layer) and the K = 0 reference (zero mean, unit scale) as out0.
import os, sys, time, importlib.util, numpy as np, flopscope as flops
from whestbench import MLP
net, spec = int(sys.argv[1]), np.load(sys.argv[2]); tag = sys.argv[3] if len(sys.argv) > 3 else "x"
M = np.asarray(spec["M"], np.float64); K = M.shape[0]
s = np.asarray(spec["s"], np.float64) if "s" in spec.files else np.ones(K)
mod_spec = importlib.util.spec_from_file_location("estf", "estimator_final_v56.py")
mod = importlib.util.module_from_spec(mod_spec); mod_spec.loader.exec_module(mod)
Wcol = np.load(f"../official/W_off{net}.npy")
mlp = MLP(width=1024, depth=16, weights=[np.ascontiguousarray(W.T).astype(np.float32) for W in Wcol])
est = mod.Estimator()


def run(m):
    est.in_mean = None if m is None else m.astype(np.float32)
    with flops.BudgetContext(flop_budget=2**45, wall_time_limit_s=3000.0, quiet=True):
        return np.asarray(est.predict(mlp, 2**41), dtype=np.float64)


t0 = time.time()
out0 = run(None)
out = np.zeros((K,) + out0.shape)
for k in range(K):
    out[k] = np.sqrt(s[k]) * run(M[k] / np.sqrt(s[k]))
    if k % 8 == 0:
        print(f"net {net}: {k + 1}/{K} at {time.time() - t0:.0f}s", flush=True)
np.savez(f"{os.environ.get('OUT', '.')}/chainat_{tag}_{net}.npz", out=out, out0=out0, M=M, s=s)
print(f"net {net}: {K} conditional chains in {time.time() - t0:.0f}s", flush=True)
