# Run the adopted system (estimator_final_v56.py, no environment needed) on one network with its per-layer dumps on,
# and save the chain's own state for theory diagnostics:
#   python chaindump.py NET
# -> $OUT/chaindump_{NET}.npz: out (16 x n predicted means), mu/var (pre-activation mean and variance per layer),
#    g4row (kappa4 diagonal per layer), C (pre-activation covariance off-diagonal per layer, float32; zero where absent)
import os, sys, importlib.util, numpy as np, flopscope as flops
from whestbench import MLP
net = int(sys.argv[1])
os.environ["V29_DUMP_LAYERS"] = ",".join(str(x) for x in range(16))
spec = importlib.util.spec_from_file_location("estf", "estimator_final_v56.py"); mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
Wcol = np.load(f"../official/W_off{net}.npy")
mlp = MLP(width=1024, depth=16, weights=[np.ascontiguousarray(W.T).astype(np.float32) for W in Wcol])
with flops.BudgetContext(flop_budget=2**45, wall_time_limit_s=3000.0, quiet=True):
    out = np.asarray(mod.Estimator().predict(mlp, 2**41), dtype=np.float64)
mt = np.load(f"../official/truth_off{net}.npz")["m"].astype(np.float64)
L, n = out.shape
mu = np.zeros((L, n)); var = np.zeros((L, n)); g4 = np.zeros((L, n)); C = np.zeros((L, n, n), np.float32)
seen = set()
for d in mod.DUMPS:
    l = d["layer"]
    if "mu" in d and l not in seen and d.get("mu") is not None:
        seen.add(l)
        mu[l] = np.asarray(d["mu"], np.float64); var[l] = np.asarray(d["var"], np.float64)
        if d.get("g4row") is not None:
            g4[l] = np.asarray(d["g4row"], np.float64)
        if d.get("C_off") is not None:
            C[l] = np.asarray(d["C_off"], np.float32)
print(f"net {net}: layers dumped {sorted(seen)}; final raw {np.mean((out[-1] - mt[-1]) ** 2):.4e}", flush=True)
np.savez(f"{os.environ.get('OUT', '.')}/chaindump_{net}.npz", out=out, mu=mu, var=var, g4row=g4, C=C)
