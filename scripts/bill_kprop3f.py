"""Billed cost (flopscope rules) and accuracy of the two-tier chain on an official network.
  /opt/wb/bin/python scripts/bill_kprop3f.py DATA NET window=2 k=128 dtype=float32 [radial=1] [collective=0]"""
import sys, os, time, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import flopscope as fs
from whest.kprop3f import kprop3f_chain
D, net = sys.argv[1], int(sys.argv[2]); o = dict(window=2, k=128, radial=1, collective=0, c4scale=1.0); dtype = "float32"
for a in sys.argv[3:]:
    k, v = a.split("=")
    if k == "dtype": dtype = v
    else: o[k] = float(v) if "." in v else int(v)
W = np.load(f"{D}/W_off{net}.npy"); mt = np.load(f"{D}/truth_off{net}.npz")["m"].astype(np.float64)
t0 = time.time()
with fs.BudgetContext(flop_budget=10 ** 14) as b:
    out = kprop3f_chain(W, o, backend="fnp", dtype=dtype)
dt = time.time() - t0; per = np.mean((out - mt) ** 2, axis=1); d = b.summary_dict()
ops = sorted(d["operations"].items(), key=lambda kv: -kv[1]["flop_cost"])[:8]
print(f"net {net} {o} {dtype}: final MSE {per[-1]:.4e}  billed C/B {b.flops_used / 2**41:.4f}  wall {dt:.0f}s (backend {d['flopscope_backend_time_s']:.0f}s, overhead {d['flopscope_overhead_time_s']:.0f}s, residual {d['residual_wall_time_s']:.1f}s)")
print("   top ops: " + ", ".join(f"{k} {v['flop_cost']/2**41:.3f}B/{v['calls']}" for k, v in ops))
print("   per layer MSE: " + " ".join(f"{p:.1e}" for p in per))
