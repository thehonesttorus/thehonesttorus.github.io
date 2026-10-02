"""Run an estimator on fetched full-split networks (N=1e9 truth); optional oracle via env H_ORACLE=<same npz>.
usage: /root/whest/bin/python run_full.py ESTIMATOR data/full_00000.npz [more.npz ...] -> prints raw (MSE vs N=1e9 truth), C/B"""
import importlib.util, os, sys, json, time
import numpy as np
import flopscope as flops
import flopscope.numpy as fnp
from whestbench import MLP, SetupContext
est_path, files = sys.argv[1], sys.argv[2:]
spec = importlib.util.spec_from_file_location('estimator', est_path); mod = importlib.util.module_from_spec(spec)
sys.path.insert(0, os.path.dirname(os.path.abspath(est_path))); spec.loader.exec_module(mod)
est = mod.Estimator(); B = 2 ** 41
with flops.BudgetContext(flop_budget=10 ** 15, quiet=True):
    est.setup(SetupContext(width=1024, depth=16, flop_budget=B, api_version='1.0', submission_dir=os.path.dirname(os.path.abspath(est_path)), seed=0))
for f in files:
    d = np.load(f); W = d['W']; T = d['truth']
    with flops.BudgetContext(flop_budget=10 ** 15, quiet=True):
        weights = [fnp.asarray(W[l]) for l in range(W.shape[0])]
    mlp = MLP(width=W.shape[1], depth=W.shape[0], weights=weights, seed=int(d['seed']) % (2 ** 63))
    ctx = flops.BudgetContext(flop_budget=B, quiet=True); t0 = time.time()
    with ctx:
        p = np.asarray(est.predict(mlp, B), dtype=np.float64); fl = ctx.flops_used
    print(json.dumps(dict(file=os.path.basename(f), raw=float(((p[-1] - T[-1]) ** 2).mean()), all_layers=float(((p - T) ** 2).mean()),
                          cb=fl / B, sec=round(time.time() - t0, 1))), flush=True)
