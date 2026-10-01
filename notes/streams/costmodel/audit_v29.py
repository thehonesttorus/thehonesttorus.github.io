"""Run 504aldo's namespace-tagged V29 (estimator_v29_ns.py, MIT) on a random He-init 1024x16 MLP
under the in-process flopscope 0.12.1 meter; dump the steady-state (second call) op ledger by
layer x family x (op_name, shapes), plus call counts and timing. Cost is data-independent."""
import sys, json, time, gc
from collections import defaultdict
sys.path.insert(0, "/tmp/claude-0/-home-user-thehonesttorus-github-io/b94ff8ab-040c-538b-a2af-dbe7aa5288dd/scratchpad/whest-p2-cumulant-k3/estimators")
import numpy as np
import flopscope as flops
import flopscope.numpy as fnp
from whestbench.domain import MLP
import estimator_v29_ns as ev

UNIT = 2.0 * 1024 ** 3
rng = np.random.default_rng(12345)
W = [(rng.standard_normal((1024, 1024)) * np.sqrt(2.0 / 1024)).astype(np.float32) for _ in range(16)]

class _Ctx:
    seed = 0
    width = 1024
    depth = 16
    flop_budget = 2 ** 41
    api_version = "1.0"
    scratch_dir = None
    submission_dir = None

def run(est, tag):
    ws = [fnp.asarray(w) for w in W]
    mlp = MLP(width=1024, depth=16, weights=ws, seed=0)
    t0 = time.time()
    with flops.BudgetContext(flop_budget=int(1e14), wall_time_limit_s=3600.0, quiet=True) as ctx:
        pred = est.predict(mlp, int(2 ** 41))
    C = float(ctx.flops_used)
    out = dict(tag=tag, units=C / UNIT, CB=C / 2 ** 41, resid=ctx.residual_wall_time_s,
               backend=ctx.flopscope_backend_time_s, overhead=ctx.flopscope_overhead_time_s,
               wall=time.time() - t0, ops=len(ctx.op_log))
    print(json.dumps(out), flush=True)
    return out, list(ctx.op_log)

est = ev.Estimator()
est.setup(_Ctx())
r1, _ = run(est, "call1")
r2, log = run(est, "call2")
lay_fam = defaultdict(lambda: defaultdict(float))
lay_fam_calls = defaultdict(lambda: defaultdict(int))
fam_op = defaultdict(lambda: defaultdict(lambda: [0.0, 0]))
for r in log:
    ns = r.namespace or ""
    segs = ns.split(".") if ns else []
    l = segs[0] if segs and segs[0].startswith("L") and segs[0][1:].isdigit() else "--"
    f = segs[-1] if segs and (len(segs) > 1 or l == "--") else "other"
    c = float(r.flop_cost)
    lay_fam[l][f] += c / UNIT
    lay_fam_calls[l][f] += 1
    key = f"{r.op_name}|{r.shapes}|{r.resolved_dtype}"
    fam_op[f][key][0] += c / UNIT
    fam_op[f][key][1] += 1
res = dict(call1=r1, call2=r2,
           layer_family={l: dict(v) for l, v in lay_fam.items()},
           layer_family_calls={l: dict(v) for l, v in lay_fam_calls.items()},
           family_ops={f: {k: v for k, v in sorted(d.items(), key=lambda kv: -kv[1][0])[:40]} for f, d in fam_op.items()})
json.dump(res, open(sys.argv[1], "w"), indent=1, default=str)
print("done", flush=True)
