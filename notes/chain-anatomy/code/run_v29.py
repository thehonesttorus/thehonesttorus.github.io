# Run the v29 estimator (or a patched copy) on official network i under the flopscope budget; print raw MSE and C/B.
import sys, os, importlib.util, numpy as np, time, flopscope as flops
from whestbench import MLP
i = int(sys.argv[1]); tag = sys.argv[2] if len(sys.argv) > 2 else "base"; src = sys.argv[3] if len(sys.argv) > 3 else "est_v29.py"
for kv in sys.argv[4:]:
    k, v = kv.split("="); os.environ[k] = v
spec = importlib.util.spec_from_file_location("estv29", src); mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
Wcol = np.load(f"../official/W_off{i}.npy"); mt = np.load(f"../official/truth_off{i}.npz")["m"].astype(np.float64)
mlp = MLP(width=1024, depth=16, weights=[np.ascontiguousarray(W.T).astype(np.float32) for W in Wcol])
est = mod.Estimator()
try:
    from whestbench import SetupContext
    est.setup(SetupContext(seed=0))
except Exception as e:
    print("setup skipped:", repr(e)[:80])
t0 = time.time()
with flops.BudgetContext(flop_budget=2**41, wall_time_limit_s=900.0, quiet=True) as ctx:
    out = np.asarray(est.predict(mlp, 2**41), dtype=np.float64)
used = getattr(ctx, "flops_used", None)
if used is None:
    try: used = ctx.summary()["flops_used"]
    except Exception: used = float("nan")
raw = np.mean((out[-1] - mt[-1])**2)
try: resid = ctx.summary()["residual_wall_time_s"]
except Exception: resid = float("nan")
print(f"net {i} {tag}: raw {raw:.4e}  C/B {used/2**41:.4f}  adjusted {raw*max(0.1, used/2**41):.4e}  wall {time.time()-t0:.0f}s  residual {resid:.3f}s", flush=True)
