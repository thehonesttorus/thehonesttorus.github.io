# Copied from notes/leaderboard-system/code/run_v29w.py (branch claude/bold-albattani-xgecxi); data path from $DATA.
# Scored-regime harness: one estimator serves a warm-up predict (another official network) and then the measured
# predict, as whestbench's worker does for every network after the first. Prints raw MSE, C/B, residual and wall time
# of the measured call.   python run_v29w.py NET TAG EST_FILE [FLAG=VAL ...]
import sys, os, importlib.util, numpy as np, time, flopscope as flops
from whestbench import MLP
i = int(sys.argv[1]); tag = sys.argv[2]; src = sys.argv[3]; DATA = os.environ.get("DATA", "../official")
for kv in sys.argv[4:]:
    k, v = kv.split("="); os.environ[k] = v
spec = importlib.util.spec_from_file_location("estv29", src); mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
def load(j):
    Wcol = np.load(f"{DATA}/W_off{j}.npy")
    return MLP(width=1024, depth=16, weights=[np.ascontiguousarray(W.T).astype(np.float32) for W in Wcol])
def _mem():
    try:
        d = dict(l.split(":", 1) for l in open("/proc/self/status") if l.startswith(("VmPeak", "VmHWM")))
        return " ".join(f"{k} {int(v.split()[0]) / 2**20:.2f}GB" for k, v in d.items())
    except Exception:
        return ""
est = mod.Estimator()
try:
    from whestbench import SetupContext
    est.setup(SetupContext(width=1024, depth=16, flop_budget=2**41, api_version="1", seed=0))
except Exception as e:
    try: est.setup(SetupContext(seed=0))
    except Exception: pass
wj = (i + 1) % 100
with flops.BudgetContext(flop_budget=2**41, wall_time_limit_s=900.0, quiet=True) as c0:
    est.predict(load(wj), 2**41)
for extra in range(int(os.environ.get("W_WARM", "1")) - 1):   # further warm-ups: measure predict #(W_WARM+1)
    with flops.BudgetContext(flop_budget=2**41, wall_time_limit_s=900.0, quiet=True):
        est.predict(load((i + 2 + extra) % 100), 2**41)
mlp = load(i); mt = np.load(f"{DATA}/truth_off{i}.npz")["m"].astype(np.float64)
t0 = time.time()
ctx = flops.BudgetContext(flop_budget=2**41, wall_time_limit_s=900.0, quiet=True)
with ctx:
    out = np.asarray(est.predict(mlp, 2**41), dtype=np.float64)
used = ctx.flops_used; raw = np.mean((out[-1] - mt[-1])**2)
if os.environ.get("SAVE_OUT"):
    np.save(os.environ["SAVE_OUT"], out.astype(np.float64))   # every layer's predicted means, for output-metric fits
print(f"net {i} {tag}: raw {raw:.4e}  C/B {used/2**41:.4f}  adjusted {raw*max(0.1, used/2**41):.4e}  wall {time.time()-t0:.1f}s  "
      f"residual {ctx.residual_wall_time_s or 0.0:.3f}s  first-call C/B {c0.flops_used/2**41:.4f} residual {c0.residual_wall_time_s or 0.0:.3f}s"
      f"  {_mem()}", flush=True)
