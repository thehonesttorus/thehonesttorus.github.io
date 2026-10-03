import sys, time, numpy as np
import flopscope as flops, flopscope.numpy as fnp
from gac_flopscope import predict
i = int(sys.argv[1]); trees = sys.argv[2] == "1"
Wcol = np.load(f"../official/W_off{i}.npy"); mt = np.load(f"../official/truth_off{i}.npz")["m"].astype(np.float64)
weights = [fnp.asarray(np.ascontiguousarray(W.T), dtype=fnp.float32) for W in Wcol]
t0 = time.time()
with flops.BudgetContext(flop_budget=2**41, wall_time_limit_s=600.0, quiet=True) as ctx:
    out = predict(weights, trees=trees)
out = np.asarray(out, dtype=np.float64)
print(f"mlp {i} trees={trees}: final MSE {np.mean((out[-1]-mt[-1])**2):.4e}  flops {ctx.flops_used:,} = {ctx.flops_used/2**41*100:.2f}% of budget  wall {time.time()-t0:.1f}s")
try: print(flops.budget_summary())
except Exception as e: print("summary:", e)
