# Free-running chain on official networks, V39_KD off/on (note XXXVI).   python run_kd.py KD NET [NET ...]
import sys, os, importlib.util, numpy as np, flopscope as flops
from whestbench import MLP
kd = sys.argv[1]; nets = [int(x) for x in sys.argv[2:]]
os.environ.update({"V29_WARM_JOIN": "1", "V29_WARM_FB": "1", "V17_R_RES": "4", "V26_STRASSEN": "6", "V26_STRASSEN_MIN": "16",
                   "V32_JOIN_SMM": "1", "V32_ROT_SMM": "1", "V32_JOIN_POST": "3", "V32_JP_C": "0.1", "V21_R_OLD": "320",
                   "V24_R_OLD2": "192", "V18_R_FB": "2", "V24_AGE_OLD2": "8", "V33_K4Q": "3", "V33_K4Q_RANK": "4",
                   "V34_OPT": "abcde", "V39_KD": kd})
spec = importlib.util.spec_from_file_location("estv29", "est_v29.py"); mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
for net in nets:
    Wcol = np.load(f"../official/W_off{net}.npy"); mt = np.load(f"../official/truth_off{net}.npz")["m"].astype(np.float64)
    mlp = MLP(width=1024, depth=16, weights=[np.ascontiguousarray(W.T).astype(np.float32) for W in Wcol])
    with flops.BudgetContext(flop_budget=2**44, wall_time_limit_s=900.0, quiet=True) as bc:
        out = np.asarray(mod.Estimator().predict(mlp, 2**41), dtype=np.float64)
    try:
        cost = bc.flops_used / 2**41
    except Exception:
        cost = float("nan")
    print(f"net {net} KD={kd}: raw {np.mean((out[-1] - mt[-1]) ** 2):.5e}  cost/B {cost:.4f}", flush=True)
