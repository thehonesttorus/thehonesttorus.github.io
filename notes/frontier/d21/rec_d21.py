"""Record the estimator's per-layer D21 (and var, D3) via its V17_DEBUG hook for one bench MLP.
usage: V17_DEBUG=1 /root/whest/bin/python rec_d21.py ESTIMATOR MLP OUT.npz"""
import importlib.util, os, sys
import numpy as np
import flopscope as flops
import flopscope.numpy as fnp
from whestbench import MLP, SetupContext
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '../../fresh-slate/bench'))
import bench
est_path, i, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
S = bench.load_set('w1024_d16')
spec = importlib.util.spec_from_file_location('estimator', est_path); mod = importlib.util.module_from_spec(spec)
sys.path.insert(0, os.path.dirname(os.path.abspath(est_path))); spec.loader.exec_module(mod)
est = mod.Estimator(); B = 2 ** 41
with flops.BudgetContext(flop_budget=10 ** 15, quiet=True):
    est.setup(SetupContext(width=1024, depth=16, flop_budget=B, api_version='1.0', submission_dir=os.path.dirname(os.path.abspath(est_path)), seed=0))
W = bench.weights(S, i)
with flops.BudgetContext(flop_budget=10 ** 15, quiet=True):
    weights = [fnp.asarray(W[l]) for l in range(W.shape[0])]
mlp = MLP(width=1024, depth=16, weights=weights, seed=int(S['estimator_seeds'][i]) % (2 ** 63))
with flops.BudgetContext(flop_budget=10 ** 15, quiet=True):
    p = np.asarray(est.predict(mlp, B), dtype=np.float64)
res = {'pred': p}
for d in mod.DEBUG:
    l = int(d['layer'])
    for k in ('D21', 'D3', 'var'):
        if d.get(k) is not None:
            res[f'{k}_{l}'] = np.asarray(d[k], dtype=np.float32)
np.savez(out, **res)
print('saved', out, sorted(k for k in res if k.startswith('D21_')))
