import sys, os, importlib.util, numpy as np, flopscope as flops
from whestbench import MLP
i = int(sys.argv[1]); scale = float(sys.argv[2]); layers = tuple(int(x) for x in sys.argv[3].split(",")) if sys.argv[3] != "none" else (); oracle = sys.argv[4] == "1"
os.environ["K3_WIN"] = "0"
spec = importlib.util.spec_from_file_location("est", "est_win.py"); mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
Wcol = np.load(f"../official/W_off{i}.npy"); mt = np.load(f"../official/truth_off{i}.npz")["m"].astype(np.float64)
mod.G4SCALE = scale; mod.G4LAYERS = layers
if oracle: mod.ORACLE_PREV = mt.astype(np.float32); mod.ORACLE_LAYERS = (15,)
mlp = MLP(width=1024, depth=16, weights=[np.ascontiguousarray(W.T).astype(np.float32) for W in Wcol])
with flops.BudgetContext(flop_budget=2**42, wall_time_limit_s=900.0, quiet=True):
    out = np.asarray(mod.Estimator().predict(mlp, 2**41), dtype=np.float64)
print(f"net {i} g4 x{scale} at layers {sys.argv[3]} oracle_mean_last={oracle}: final MSE {np.mean((out[-1]-mt[-1])**2):.3e}", flush=True)
