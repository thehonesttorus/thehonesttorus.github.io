# Cost profile of the chain by source line: every billed op is attributed to (a) the innermost est file line and
# (b) the call site inside _predict_core (helpers and the Strassen class roll up to it).
#   python prof_chain.py NET EST_FILE "FLAGS"
import sys, os, collections, linecache, importlib.util, numpy as np
net, est = int(sys.argv[1]), sys.argv[2]
for kv in sys.argv[3].split():
    k, v = kv.split("="); os.environ[k] = v
import flopscope as flops
import flopscope._budget as B
from whestbench import MLP
cls = [c for c in vars(B).values() if isinstance(c, type) and hasattr(c, "_append_op_record")][0]
estpath = os.path.abspath(est)
inner, outer, byop = collections.Counter(), collections.Counter(), collections.Counter()
cnt_in, shapes_in = collections.Counter(), {}
orig = cls._append_op_record
def rec(self, record):
    f = sys._getframe(1); il = ol = None; fn_outer = None
    while f is not None:
        if f.f_code.co_filename == estpath:
            if il is None: il = (f.f_lineno, f.f_code.co_name)
            if f.f_code.co_name == "_predict_core": ol = f.f_lineno; break
        f = f.f_back
    c = record.flop_cost
    inner[il] += c; outer[ol] += c; byop[(record.op_name, record.resolved_dtype)] += c
    cnt_in[il] += 1; shapes_in.setdefault(il, (record.op_name, record.shapes))
    return orig(self, record)
cls._append_op_record = rec
spec = importlib.util.spec_from_file_location("estmod", est); mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
Wcol = np.load(f"../official/W_off{net}.npy"); mt = np.load(f"../official/truth_off{net}.npz")["m"].astype(np.float64)
mlp = MLP(width=1024, depth=16, weights=[np.ascontiguousarray(W.T).astype(np.float32) for W in Wcol])
with flops.BudgetContext(flop_budget=2**44, wall_time_limit_s=900.0, quiet=True) as bc:
    out = np.asarray(mod.Estimator().predict(mlp, 2**41), dtype=np.float64)
    total = bc.flops_used
U = 2 * 1024 ** 3
print(f"net {net}: raw {np.mean((out[-1] - mt[-1])**2):.4e}  total {total:.4e} FLOPs = {total / 2**41:.4f} B = {total / U:.1f} units (2n^3)")
def src(line): return linecache.getline(estpath, line).strip()[:110] if line else "(outside est)"
print("\n== by call site in _predict_core (units, % of total)")
for l, c in outer.most_common(45):
    print(f"{c / U:8.2f} {100 * c / total:5.1f}%  L{l}: {src(l)}")
print("\n== by innermost line (units, count, first op)")
for il, c in inner.most_common(60):
    l, fn = il if il else (None, None)
    print(f"{c / U:8.2f} {100 * c / total:5.1f}% n={cnt_in[il]:6d} {fn}:L{l} {shapes_in[il][0]} {str(shapes_in[il][1])[:60]} | {src(l)}")
print("\n== by op and dtype")
for (op, dt), c in byop.most_common(30):
    print(f"{c / U:8.2f} {100 * c / total:5.1f}%  {op} [{dt}]")
