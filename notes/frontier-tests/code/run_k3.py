# Run the open-source K=3 estimator on official networks with flopscope metering.
import sys, time, importlib.util, numpy as np
import flopscope as flops, flopscope.numpy as fnp
spec = importlib.util.spec_from_file_location("est", sys.argv[1]); mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
from whestbench import MLP
est = mod.Estimator()
try:
    from whestbench import SetupContext
    est.setup(SetupContext(seed=0)) if hasattr(est, "setup") else None
except Exception as e: print("setup:", e)
for i in [int(x) for x in sys.argv[2].split(",")]:
    Wcol = np.load(f"../official/W_off{i}.npy"); mt = np.load(f"../official/truth_off{i}.npz")["m"].astype(np.float64)
    weights = [np.ascontiguousarray(W.T).astype(np.float32) for W in Wcol]
    mlp = MLP(width=1024, depth=16, weights=weights)
    t0 = time.time()
    with flops.BudgetContext(flop_budget=2**41, wall_time_limit_s=600.0, quiet=True) as ctx:
        out = est.predict(mlp, 2**41)
    out = np.asarray(out, dtype=np.float64)
    print(f"mlp {i}: final MSE {np.mean((out[-1]-mt[-1])**2):.4e}  C/B {ctx.flops_used/2**41:.4f}  wall {time.time()-t0:.1f}s", flush=True)
