import sys, os, importlib.util, numpy as np, flopscope as flops
from whestbench import MLP
sys.path.insert(0, "../num12")
from gac import gac
i = int(sys.argv[1]); layers = tuple(int(x) for x in sys.argv[2].split(",")); gscale = float(sys.argv[3])
os.environ["K3_WIN"] = "0"
spec = importlib.util.spec_from_file_location("est", "est_win.py"); mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
Wcol = np.load(f"../official/W_off{i}.npy"); mt = np.load(f"../official/truth_off{i}.npz")["m"].astype(np.float64)
gam = np.array([d["gamma"] for d in gac(list(Wcol.astype(np.float64)))])
mod.GAIN_K4 = np.maximum(gam*gscale, 1e-6); mod.G4LAYERS = layers; mod.G4SCALE = 1.0
mlp = MLP(width=1024, depth=16, weights=[np.ascontiguousarray(W.T).astype(np.float32) for W in Wcol])
with flops.BudgetContext(flop_budget=2**42, wall_time_limit_s=900.0, quiet=True):
    out = np.asarray(mod.Estimator().predict(mlp, 2**41), dtype=np.float64)
print(f"net {i} gain-k4 (gamma x{gscale}) at layers {sys.argv[2]}: final MSE {np.mean((out[-1]-mt[-1])**2):.3e}", flush=True)
