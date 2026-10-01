"""Stage P quick runner: a flopscope estimator (whest Estimator class in estimator.py) run in-process on a bench set,
weights regenerated from the seeds (no parquet needed), metered under BudgetContext(2^41) exactly as the local runner.
Reports per-MLP final MSE, raw (minus truth noise), C/B, residual wall time, wall; adjusted = raw * max(0.1, C/B).
It does NOT enforce the grader's process isolation or 8 GB limit: use run_p.py (`whest run --runner subprocess`) for the
final check.

  python run_p_inproc.py --estimator path/estimator.py --set w1024_d16 [--mlps 0,1] [--json out.json]
"""
import argparse, importlib.util, json, os, sys, time

ap = argparse.ArgumentParser()
ap.add_argument("--estimator", required=True)
ap.add_argument("--set", required=True)
ap.add_argument("--mlps", default="")
ap.add_argument("--json")
a = ap.parse_args()

import numpy as np
import flopscope as flops
import flopscope.numpy as fnp
from whestbench import MLP, SetupContext

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bench  # noqa: E402

S = bench.load_set(a.set)
idx = [int(x) for x in a.mlps.split(",")] if a.mlps else list(range(len(S["seeds"])))
spec = importlib.util.spec_from_file_location("estimator", a.estimator)
mod = importlib.util.module_from_spec(spec)
sys.path.insert(0, os.path.dirname(os.path.abspath(a.estimator)))
spec.loader.exec_module(mod)
cls = getattr(mod, "Estimator")
est = cls()
B = 2 ** 41
t0 = time.time()
with flops.BudgetContext(flop_budget=10 ** 15, quiet=True):
    est.setup(SetupContext(width=S["width"], depth=S["depth"], flop_budget=B, api_version="1.0",
                           submission_dir=os.path.dirname(os.path.abspath(a.estimator)), seed=0))
print(f"setup {time.time() - t0:.2f}s")
rows = []
for i in idx:
    W = bench.weights(S, i)
    with flops.BudgetContext(flop_budget=10 ** 15, quiet=True):
        weights = [fnp.asarray(W[l]) for l in range(W.shape[0])]
    mlp = MLP(width=S["width"], depth=S["depth"], weights=weights, seed=int(S["estimator_seeds"][i]) % (2 ** 63))
    ctx = flops.BudgetContext(flop_budget=B, quiet=True)
    rec = dict(mlp=i)
    t1 = time.time()
    try:
        with ctx:
            pred = est.predict(mlp, B)
            fl = ctx.flops_used
        p = np.asarray(pred, dtype=np.float64)
        T = S["means"][i]
        mse = float(np.mean((p[-1] - T[-1]) ** 2))
        raw = mse - float(S["noise"][i])
        rec.update(ok=bool(np.isfinite(p).all()), cb=fl / B, mse=mse, raw=raw, adjusted=raw * max(0.1, fl / B),
                   all_layers=float(np.mean((p - T) ** 2)))
    except Exception as e:  # noqa: BLE001
        rec.update(ok=False, error=f"{type(e).__name__}: {e}"[:300])
    rec.update(wall=time.time() - t1, residual=ctx.residual_wall_time_s)
    rows.append(rec)
    print(json.dumps(rec), flush=True)
ok = [r for r in rows if r.get("ok")]
if ok:
    print(f"SUMMARY {a.set}: {len(ok)}/{len(rows)} ok, raw {np.mean([r['raw'] for r in ok]):.4e}, "
          f"C/B max {max(r['cb'] for r in ok):.4f}, adjusted {np.mean([r['adjusted'] for r in ok]):.4e}, "
          f"residual max {max(r['residual'] for r in rows):.3f}s, wall max {max(r['wall'] for r in rows):.1f}s")
if a.json:
    json.dump(dict(set=a.set, estimator=a.estimator, per_mlp=rows), open(a.json, "w"), indent=1)
