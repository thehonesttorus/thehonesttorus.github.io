# Late-layer oracle ceiling (note XLII, N3 referee's gate): the adopted system with some of its readouts replaced by
# Monte Carlo truth at chosen layers only.
#   V37_ORACLE=D3,G4 V37_ORACLE_LAYERS=12,13,14,15 V37_ORACLE_FILE=mc2_off0_full.npz python oracle_late.py NET TAG
# Prints the per-layer MSE against the dataset truth and saves $OUT/oracle_{TAG}_{NET}.npy (16 x n means).
import os, sys, importlib.util, numpy as np, flopscope as flops
from whestbench import MLP
net, tag = int(sys.argv[1]), sys.argv[2]
spec = importlib.util.spec_from_file_location("estf", "estimator_final_v56.py")
mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
Wcol = np.load(f"../official/W_off{net}.npy")
mlp = MLP(width=1024, depth=16, weights=[np.ascontiguousarray(W.T).astype(np.float32) for W in Wcol])
with flops.BudgetContext(flop_budget=2**45, wall_time_limit_s=3000.0, quiet=True):
    out = np.asarray(mod.Estimator().predict(mlp, 2**41), dtype=np.float64)
mt = np.load(f"../official/truth_off{net}.npz")["m"].astype(np.float64)
mse = np.mean((out - mt) ** 2, axis=1)
print(f"net {net} {tag} oracle={os.environ.get('V37_ORACLE', '')} layers={os.environ.get('V37_ORACLE_LAYERS', 'all')} "
      f"file={os.environ.get('V37_ORACLE_FILE', '')}: final raw {mse[-1]:.5e} | " + " ".join(f"{l}:{mse[l]:.3e}" for l in (9, 12, 13, 14, 15)),
      flush=True)
np.save(f"{os.environ.get('OUT', '.')}/oracle_{tag}_{net}.npy", out)
