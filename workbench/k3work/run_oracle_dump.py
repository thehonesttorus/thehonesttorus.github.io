# Chain with the TRUE post-activation mean injected entering every layer; dump the chain's per-layer state.
import sys, os, importlib.util, numpy as np, pickle, flopscope as flops
from whestbench import MLP
i = int(sys.argv[1])
os.environ["K3_WIN"] = "0"; os.environ["K3_DUMP"] = "1"
spec = importlib.util.spec_from_file_location("est", "est_win.py"); mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
Wcol = np.load(f"../official/W_off{i}.npy"); mt = np.load(f"../official/truth_off{i}.npz")["m"].astype(np.float64)
mod.ORACLE_PREV = mt.astype(np.float32); mod.ORACLE_LAYERS = tuple(range(1, 16))
mlp = MLP(width=1024, depth=16, weights=[np.ascontiguousarray(W.T).astype(np.float32) for W in Wcol])
with flops.BudgetContext(flop_budget=2**42, wall_time_limit_s=900.0, quiet=True):
    out = np.asarray(mod.Estimator().predict(mlp, 2**41), dtype=np.float64)
print("local MSE per layer:", " ".join(f"{np.mean((out[l]-mt[l])**2):.2e}" for l in range(16)), flush=True)
np.save(f"pred_oracleall_off{i}.npy", out)
pickle.dump(mod.DUMP, open(f"dump_oracleall_off{i}.pkl", "wb"))
