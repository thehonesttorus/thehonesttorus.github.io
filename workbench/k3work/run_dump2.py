# Plain chain run (no oracle) with the per-layer dump and the last-layer g4row dump.
import sys, os, importlib.util, numpy as np, pickle, flopscope as flops
from whestbench import MLP
i = int(sys.argv[1])
os.environ["K3_WIN"] = "0"; os.environ["K3_DUMP"] = "1"; os.environ["K3_DUMP2"] = "1"
spec = importlib.util.spec_from_file_location("est", "est_win.py"); mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
Wcol = np.load(f"../official/W_off{i}.npy"); mt = np.load(f"../official/truth_off{i}.npz")["m"].astype(np.float64)
mlp = MLP(width=1024, depth=16, weights=[np.ascontiguousarray(W.T).astype(np.float32) for W in Wcol])
with flops.BudgetContext(flop_budget=2**42, wall_time_limit_s=900.0, quiet=True) as ctx:
    out = np.asarray(mod.Estimator().predict(mlp, 2**41), dtype=np.float64)
print("final MSE", np.mean((out[-1]-mt[-1])**2), "C/B", ctx.flops_used/2**41 if hasattr(ctx,'flops_used') else None, flush=True)
pickle.dump(dict(DUMP=mod.DUMP, DUMP2=mod.DUMP2), open(f"dump2_off{i}.pkl", "wb"))
