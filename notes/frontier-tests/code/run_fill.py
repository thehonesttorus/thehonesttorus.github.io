# Collar + gain fill: youngest WIN sources exact, dropped sources replaced by the scale-mixture third cumulant
#   kappa3_old(i,j,k) = (g_old/2)(mu_i S_jk + mu_j S_ik + mu_k S_ij),  g_old(l) = sum of weights-only gain injections
#   born in the dropped layers (critical: no decay).
import sys, os, time, importlib.util, numpy as np
import flopscope as flops
from whestbench import MLP
sys.path.insert(0, "../num12")
from gac import gac
i, win, fill = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3] == "1"
os.environ["K3_WIN"] = str(win); os.environ["K3_DUMP"] = "0"
spec = importlib.util.spec_from_file_location("est", "est_win.py"); mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
Wcol = np.load(f"../official/W_off{i}.npy"); mt = np.load(f"../official/truth_off{i}.npz")["m"].astype(np.float64)
gam = [d["gamma"] for d in gac(list(Wcol.astype(np.float64)))]; f = np.diff([0.0] + gam)   # f[b]: injection seen at pre-activation layer b
gold = np.zeros(16)
for li in range(16):
    nb = li - win          # sources born at layers 0..li-win-1 are dropped -> injections at pre-activation layers 1..li-win
    if win > 0 and nb >= 1: gold[li] = f[1:nb+1].sum()
mod.GOLD = gold*float(os.environ.get("FILL_SCALE", "1")) if fill else None
mlp = MLP(width=1024, depth=16, weights=[np.ascontiguousarray(W.T).astype(np.float32) for W in Wcol])
with flops.BudgetContext(flop_budget=2**42, wall_time_limit_s=900.0, quiet=True) as ctx:
    out = np.asarray(mod.Estimator().predict(mlp, 2**41), dtype=np.float64)
print(f"mlp {i} WIN={win} fill={int(fill)}: final MSE {np.mean((out[-1]-mt[-1])**2):.4e}   (g_old at last layer {gold[-1]:.4f})", flush=True)
