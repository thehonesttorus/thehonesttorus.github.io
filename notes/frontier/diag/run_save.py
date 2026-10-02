"""Run a flopscope estimator in-process on bench MLPs and save the full (L, n) predictions next to the truth.
usage: /root/whest/bin/python run_save.py ESTIMATOR MLPS OUT.npz"""
import importlib.util, os, sys, time
import numpy as np
import flopscope as flops
import flopscope.numpy as fnp
from whestbench import MLP, SetupContext
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '../../fresh-slate/bench'))
import bench
est_path, mlps, out = sys.argv[1], [int(a) for a in sys.argv[2].split(',')], sys.argv[3]
S = bench.load_set('w1024_d16')
spec = importlib.util.spec_from_file_location('estimator', est_path); mod = importlib.util.module_from_spec(spec)
sys.path.insert(0, os.path.dirname(os.path.abspath(est_path))); spec.loader.exec_module(mod)
est = mod.Estimator(); B = 2 ** 41
with flops.BudgetContext(flop_budget=10 ** 15, quiet=True):
    est.setup(SetupContext(width=S['width'], depth=S['depth'], flop_budget=B, api_version='1.0',
                           submission_dir=os.path.dirname(os.path.abspath(est_path)), seed=0))
P, T, NZ, CB = [], [], [], []
for i in mlps:
    W = bench.weights(S, i)
    with flops.BudgetContext(flop_budget=10 ** 15, quiet=True):
        weights = [fnp.asarray(W[l]) for l in range(W.shape[0])]
    mlp = MLP(width=S['width'], depth=S['depth'], weights=weights, seed=int(S['estimator_seeds'][i]) % (2 ** 63))
    ctx = flops.BudgetContext(flop_budget=B, quiet=True)
    t0 = time.time()
    with ctx:
        p = np.asarray(est.predict(mlp, B), dtype=np.float64)
        fl = ctx.flops_used
    P.append(p); T.append(S['means'][i]); NZ.append(S['noise'][i]); CB.append(fl / B)
    print(i, 'raw %.4e' % (((p[-1] - S['means'][i][-1]) ** 2).mean() - S['noise'][i]), 'cb %.4f' % (fl / B), '%.0fs' % (time.time() - t0), flush=True)
np.savez(out, pred=np.array(P), truth=np.array(T), noise=np.array(NZ), cb=np.array(CB), mlps=np.array(mlps))
