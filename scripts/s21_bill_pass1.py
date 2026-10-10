"""Bill and check the flopscope pass 1 (whest/s21_fnp.py) on an official network, exactly as the grader would meter it.
/opt/wb/bin/python scripts/s21_bill_pass1.py DATA NET [KCG_FILE]"""
import sys, os, time, collections
import numpy as np
import flopscope as flops, flopscope.numpy as fnp
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from whest.s21_fnp import pass1
D, net = sys.argv[1], int(sys.argv[2]); B = 2 ** 41
weights = [fnp.asarray(np.ascontiguousarray(W.T).astype(np.float32)) for W in np.load(f"{D}/W_off{net}.npy")]
truth = np.load(f"{D}/truth_off{net}.npz")["m"].astype(np.float64)
with flops.BudgetContext(flop_budget=B, quiet=True):            # warm-up (path caches), as the grader's worker would
    pass1(weights)
t0 = time.time()
with flops.BudgetContext(flop_budget=B, quiet=True) as ctx:
    out, gaps, _ = pass1(weights)
    out = fnp.asarray(out, dtype=fnp.float32)
wall = time.time() - t0
ops = collections.Counter(); calls = collections.Counter()
for r in ctx.op_log: ops[(r.op_name, str(r.resolved_dtype))] += r.flop_cost; calls[r.op_name] += 1
y = np.asarray(out, dtype=np.float64)
print(f"net {net}: billed {ctx.flops_used:.4e} = {ctx.flops_used / 2**30:.1f} n^3 = {ctx.flops_used / B:.4f} B | counted calls {len(ctx.op_log)} |"
      f" wall {wall:.2f}s residual {ctx.residual_wall_time_s:.3f}s backend {ctx.flopscope_backend_time_s:.2f}s overhead {ctx.flopscope_overhead_time_s:.2f}s")
print("  top ops: " + "; ".join(f"{k[0]}[{k[1]}] {v / 2**30:.2f} n^3" for k, v in ops.most_common(8)))
print(f"  final-layer RMS vs truth {np.sqrt(np.mean((y[-1] - truth[-1]) ** 2)):.3e}")
if len(sys.argv) > 3:
    ref = np.load(sys.argv[3])["m"]
    print(f"  f32 flopscope vs f64 numpy closure: final-layer RMS diff {np.sqrt(np.mean((y[-1] - ref[-1]) ** 2)):.2e}, max over layers "
          f"{max(np.sqrt(np.mean((y[l] - ref[l]) ** 2)) for l in range(len(y))):.2e}")
print("  gaps g_l: " + " ".join(f"{float(g):.3f}" for g in gaps))
