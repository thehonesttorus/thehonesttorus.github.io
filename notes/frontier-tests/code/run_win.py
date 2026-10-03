# Full chain vs youngest-WIN-sources chain on official networks; dump per-layer D3/D21/mu/C_pre.
import sys, os, time, importlib.util, numpy as np, pickle
import flopscope as flops, flopscope.numpy as fnp
from whestbench import MLP
i = int(sys.argv[1]); win = sys.argv[2]
os.environ["K3_WIN"] = win; os.environ["K3_DUMP"] = "1"
spec = importlib.util.spec_from_file_location("est", "est_win.py"); mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
Wcol = np.load(f"../official/W_off{i}.npy"); mt = np.load(f"../official/truth_off{i}.npz")["m"].astype(np.float64)
mlp = MLP(width=1024, depth=16, weights=[np.ascontiguousarray(W.T).astype(np.float32) for W in Wcol])
t0 = time.time()
with flops.BudgetContext(flop_budget=2**42, wall_time_limit_s=900.0, quiet=True) as ctx:
    out = np.asarray(mod.Estimator().predict(mlp, 2**41), dtype=np.float64)
print(f"mlp {i} WIN={win}: final MSE {np.mean((out[-1]-mt[-1])**2):.4e}  C/B {ctx.flops_used/2**41:.3f}  {time.time()-t0:.0f}s", flush=True)
np.save(f"pred_off{i}_win{win}.npy", out)
pickle.dump([{k: v for k, v in d.items()} for d in mod.DUMP], open(f"dump_off{i}_win{win}.pkl", "wb"))
