# Grader reproduction: whestbench's own SubprocessRunner (the evaluator's worker: RLIMIT_AS = memory_limit_mb,
# setup() outside the meter, one estimator instance serving every network in order).
#   python run_worker.py NETS EST_FILE MEM_MB [FLAG=VAL ...]      NETS = "0-5" or "3,7,9"
# Per network prints status, raw MSE, C/B, residual, wall; the worker's peak memory comes from /proc of the child.
import sys, os, time, threading
import numpy as np
nets_s, est, mem = sys.argv[1], sys.argv[2], int(sys.argv[3])
for kv in sys.argv[4:]:
    k, v = kv.split("=", 1); os.environ[k] = v
from pathlib import Path
from whestbench.runner import SubprocessRunner, EstimatorEntrypoint, ResourceLimits
from whestbench.sdk import SetupContext
from whestbench import MLP
if "-" in nets_s:
    a, b = map(int, nets_s.split("-")); nets = list(range(a, b + 1))
else:
    nets = [int(x) for x in nets_s.split(",")]
lim = ResourceLimits(setup_timeout_s=5.0, predict_timeout_s=600.0, memory_limit_mb=mem, flop_budget=2 ** 41,
                     wall_time_limit_s=120.0, residual_wall_time_limit_s=0.4)
r = SubprocessRunner()
peak = {"hwm": 0.0, "vm": 0.0}
def watch():
    while True:
        p = getattr(r, "_process", None)
        if p is not None and p.poll() is None:
            try:
                for l in open(f"/proc/{p.pid}/status"):
                    if l.startswith("VmHWM"): peak["hwm"] = max(peak["hwm"], int(l.split()[1]) / 2 ** 20)
                    if l.startswith("VmPeak"): peak["vm"] = max(peak["vm"], int(l.split()[1]) / 2 ** 20)
            except Exception:
                pass
        time.sleep(0.2)
threading.Thread(target=watch, daemon=True).start()
t0 = time.time()
try:
    r.start(EstimatorEntrypoint(file_path=Path(est).resolve(), class_name="Estimator"),
            SetupContext(width=1024, depth=16, flop_budget=2 ** 41, api_version="1.0", seed=0), lim)
    print(f"setup ok {time.time() - t0:.2f}s", flush=True)
except Exception as e:
    print(f"setup FAILED {type(e).__name__}: {str(e)[:400]}", flush=True); sys.exit(1)
for i in nets:
    Wcol = np.load(f"../official/W_off{i}.npy")
    mlp = MLP(width=1024, depth=16, weights=[np.ascontiguousarray(W.T).astype(np.float32) for W in Wcol])
    t1 = time.time()
    try:
        out = np.asarray(r.predict(mlp, 2 ** 41), dtype=np.float64)
        st = r.last_predict_stats()
        mt = np.load(f"../official/truth_off{i}.npz")["m"].astype(np.float64)
        raw = np.mean((out[-1] - mt[-1]) ** 2)
        print(f"net {i}: ok raw {raw:.4e} C/B {st.flops_used / 2 ** 41:.4f} residual {st.residual_wall_time_s:.3f}s "
              f"wall {st.wall_time_s:.1f}s | worker VmHWM {peak['hwm']:.2f}GB VmPeak {peak['vm']:.2f}GB", flush=True)
    except Exception as e:
        print(f"net {i}: FAILED after {time.time() - t1:.1f}s {type(e).__name__}: {str(e)[:600]} | worker VmHWM {peak['hwm']:.2f}GB "
              f"VmPeak {peak['vm']:.2f}GB", flush=True)
r.close()
