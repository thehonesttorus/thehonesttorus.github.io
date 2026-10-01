"""Profile an estimator under the in-process flopscope 0.12.1 meter: units, calls and residual
per family (namespace leaf) and per layer, steady state (the last of --calls predicts).
Residual is attributed to ops as the wall gap since the previous op ended
(start_i - start_{i-1} - backend_{i-1} - overhead_{i-1}), i.e. the user Python before op i.
Usage: prof.py EST.py OUT.json [--mlp he|dev:IDX] [--calls 2] [--save-pred PRED.npy] [--threads 2]
"""
import argparse, gc, importlib.util, json, os, sys, time
from collections import defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("est"); ap.add_argument("out")
ap.add_argument("--mlp", default="he")
ap.add_argument("--calls", type=int, default=2)
ap.add_argument("--save-pred", default=None)
ap.add_argument("--dataset", default="/root/work/dev6")
ap.add_argument("--keep-log", action="store_true")
a = ap.parse_args()

import numpy as np
import flopscope as flops
import flopscope.numpy as fnp
from whestbench.domain import MLP

UNIT = 2.0 * 1024 ** 3
if a.mlp == "he":
    rng = np.random.default_rng(12345)
    Wn = [(rng.standard_normal((1024, 1024)) * np.sqrt(2.0 / 1024)).astype(np.float32) for _ in range(16)]
    seed = 0
else:
    import glob, pyarrow.parquet as pq
    idx = int(a.mlp.split(":")[1])
    f = sorted(glob.glob(os.path.join(a.dataset, "data", "*.parquet")))[0]
    r = pq.read_table(f).slice(idx, 1).to_pylist()[0]
    w = np.asarray(r["weights"], dtype=np.float32)
    Wn = [w[l] for l in range(w.shape[0])]
    seed = int(r["mlp_seed"]) % (2 ** 63)

spec = importlib.util.spec_from_file_location("est_mod", a.est)
mod = importlib.util.module_from_spec(spec)
sys.path.insert(0, os.path.dirname(os.path.abspath(a.est)))
spec.loader.exec_module(mod)


class _Ctx:
    seed = 0; width = len(Wn[0]); depth = len(Wn); flop_budget = 2 ** 41
    api_version = "1.0"; scratch_dir = None; submission_dir = os.path.dirname(os.path.abspath(a.est))


est = mod.Estimator()
with flops.BudgetContext(flop_budget=10 ** 15, quiet=True):
    est.setup(_Ctx())
res = {"est": a.est, "mlp": a.mlp, "calls": []}
for c in range(a.calls):
    with flops.BudgetContext(flop_budget=10 ** 15, quiet=True):
        ws = [fnp.asarray(x) for x in Wn]
    mlp = MLP(width=len(Wn[0]), depth=len(Wn), weights=ws, seed=seed)
    t0 = time.time()
    ctx = flops.BudgetContext(flop_budget=2 ** 41, wall_time_limit_s=3600.0, quiet=True)
    with ctx:
        pred = est.predict(mlp, 2 ** 41)
    rec = dict(units=ctx.flops_used / UNIT, CB=ctx.flops_used / 2 ** 41, resid=ctx.residual_wall_time_s,
               backend=ctx.flopscope_backend_time_s, overhead=ctx.flopscope_overhead_time_s,
               wall=time.time() - t0, ops=len(ctx.op_log))
    print(json.dumps(rec), flush=True)
    res["calls"].append(rec)
    log = list(ctx.op_log)
    pred_np = np.asarray(pred)
    if a.save_pred:
        np.save(a.save_pred.replace(".npy", f"_c{c}.npy"), pred_np)
    del ctx
    gc.collect()

fam = defaultdict(lambda: [0.0, 0, 0.0])
lay = defaultdict(lambda: [0.0, 0, 0.0])
opk = defaultdict(lambda: [0.0, 0, 0.0])
prev_end = 0.0
gsum = 0.0
for rr in log:
    ns = rr.namespace or ""
    segs = ns.split(".") if ns else []
    l = segs[0] if segs and segs[0].startswith("L") and segs[0][1:].isdigit() else "--"
    f = segs[-1] if segs and (len(segs) > 1 or l == "--") else "other"
    st = rr.flopscope_context_start_offset_s or 0.0
    gap = max(0.0, st - prev_end)
    prev_end = st + (rr.flopscope_backend_duration_s or 0.0) + (rr.flopscope_overhead_duration_s or 0.0)
    gsum += gap
    u = float(rr.flop_cost) / UNIT
    for d, key in ((fam, f), (lay, l), (opk, f"{f}|{rr.op_name}")):
        d[key][0] += u; d[key][1] += 1; d[key][2] += gap
res["gap_sum"] = gsum
res["family"] = {k: dict(units=v[0], calls=v[1], gap_ms=1e3 * v[2]) for k, v in sorted(fam.items(), key=lambda kv: -kv[1][0])}
res["layer"] = {k: dict(units=v[0], calls=v[1], gap_ms=1e3 * v[2]) for k, v in sorted(lay.items())}
res["family_op"] = {k: dict(units=v[0], calls=v[1], gap_ms=1e3 * v[2]) for k, v in sorted(opk.items(), key=lambda kv: -kv[1][2])[:80]}
json.dump(res, open(a.out, "w"), indent=1)
print(f"gap_sum {gsum:.3f}s  resid {res['calls'][-1]['resid']:.3f}s")
for k, v in list(res["family"].items()):
    print(f"{k:18s} {v['units']:8.2f} u {v['calls']:6d} calls {v['gap_ms']:8.1f} ms")
