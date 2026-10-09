# Cold-screen run of the adopted system (estimator_final_v56.py) under whatever research switches the environment sets.
#   [V58_K31SC=1 ...] python screen56.py NET TAG
# Prints the final-layer raw MSE against the dataset truth and the flopscope cost; with SAVE_OUT, saves all 16 layers'
# predicted means.
import os, sys, importlib.util, numpy as np, flopscope as flops
from whestbench import MLP
net, tag = int(sys.argv[1]), sys.argv[2]
spec = importlib.util.spec_from_file_location("estf", "estimator_final_v56.py"); mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
OFF = os.environ.get("OFFDIR", "../official")
Wcol = np.load(f"{OFF}/W_off{net}.npy"); mt = np.load(f"{OFF}/truth_off{net}.npz")["m"].astype(np.float64)
mlp = MLP(width=1024, depth=16, weights=[np.ascontiguousarray(W.T).astype(np.float32) for W in Wcol])
with flops.BudgetContext(flop_budget=2**44, wall_time_limit_s=900.0, quiet=True) as _bc:
    out = np.asarray(mod.Estimator().predict(mlp, 2**41), dtype=np.float64)
_F = float(_bc.flops_used() if callable(_bc.flops_used) else _bc.flops_used)
_R = float(_bc.residual_wall_time_s() if callable(_bc.residual_wall_time_s) else _bc.residual_wall_time_s)
if os.environ.get("SAVE_OUT"):
    np.save(os.environ["SAVE_OUT"], out)
print(f"net {net} {tag}: raw {np.mean((out[-1] - mt[-1]) ** 2):.5e} | F {_F:.5e} ({_F / 2**41:.5f} B) residual wall "
      f"{_R:.3f} s", flush=True)
