# Oracle: replace the chain's propagated post-activation mean entering layer(s) L by the TRUE mean (from 1e9-sample truth).
import sys, os, importlib.util, numpy as np, flopscope as flops
from whestbench import MLP
i = int(sys.argv[1]); layers = tuple(int(x) for x in sys.argv[2].split(",")) if sys.argv[2] != "none" else ()
os.environ["K3_WIN"] = "0"; os.environ["K3_DUMP"] = "0"
spec = importlib.util.spec_from_file_location("est", "est_win.py"); mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
Wcol = np.load(f"../official/W_off{i}.npy"); mt = np.load(f"../official/truth_off{i}.npz")["m"].astype(np.float64)
mod.ORACLE_PREV = mt.astype(np.float32); mod.ORACLE_LAYERS = layers
mlp = MLP(width=1024, depth=16, weights=[np.ascontiguousarray(W.T).astype(np.float32) for W in Wcol])
with flops.BudgetContext(flop_budget=2**42, wall_time_limit_s=900.0, quiet=True):
    out = np.asarray(mod.Estimator().predict(mlp, 2**41), dtype=np.float64)
print(f"net {i}, true mean injected entering layers {layers}: final MSE {np.mean((out[-1]-mt[-1])**2):.3e}; layer MSEs " + " ".join(f"{np.mean((out[l]-mt[l])**2):.1e}" for l in [3, 7, 11, 14]), flush=True)
