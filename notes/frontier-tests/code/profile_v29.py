import sys, time, importlib.util, numpy as np, json
import flopscope as flops
spec = importlib.util.spec_from_file_location("est", sys.argv[1]); mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
from whestbench import MLP, SetupContext
est = mod.Estimator()
try: est.setup(SetupContext(seed=0))
except Exception as e: print("setup:", e)
i = int(sys.argv[2])
Wcol = np.load(f"../official/W_off{i}.npy"); mt = np.load(f"../official/truth_off{i}.npz")["m"].astype(np.float64)
mlp = MLP(width=1024, depth=16, weights=[np.ascontiguousarray(W.T).astype(np.float32) for W in Wcol])
t0 = time.time()
with flops.BudgetContext(flop_budget=2**41, wall_time_limit_s=900.0, quiet=True) as ctx:
    out = est.predict(mlp, 2**41)
out = np.asarray(out, dtype=np.float64)
print(f"net {i}: final MSE {np.mean((out[-1]-mt[-1])**2):.4e}  C/B {ctx.flops_used/2**41:.4f}  wall {time.time()-t0:.1f}s", flush=True)
d = ctx.summary_dict()
print(json.dumps({k: v for k, v in d.items() if not isinstance(v, (list, dict))}, indent=1, default=str)[:2000])
ops = d.get("ops") or d.get("operations") or d.get("by_op") or {}
print("keys:", list(d.keys()))
print(ctx.summary()[:6000])
