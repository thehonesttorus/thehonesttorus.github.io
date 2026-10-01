import sys, json, numpy as np
sys.path.insert(0, '../../bench'); sys.path.insert(0, '.')
import eval_q, bethe
modes = {'tree': bethe.estimate_tree, 'gauss': bethe.estimate_gauss,
         'fact0': lambda W: bethe.estimate_fact(W, old=0), 'fact1': lambda W: bethe.estimate_fact(W, old=1),
         'full': bethe.estimate_full,
         'edge0': lambda W: bethe.estimate_edge(W, old=0), 'edge1': lambda W: bethe.estimate_edge(W, old=1),
         'edgeG': lambda W: bethe.estimate_edge(W, old=0, hubs=False)}
units = {'tree': 1, 'gauss': 32, 'fact0': 105, 'fact1': 195, 'full': 195, 'edge0': 105, 'edge1': 195, 'edgeG': 40}
sel = sys.argv[1].split(','); sets = sys.argv[2].split(',')
out = {}
for k in sel:
    print('=====', k, flush=True)
    r = eval_q.evaluate(lambda W: modes[k](np.asarray(W, np.float64)), sets, units_1024=units[k], verbose=True)
    out[k] = dict(fit=r['fit'], sets=[{kk: v for kk, v in s.items() if kk != 'per_mlp'} for s in r['sets']])
    json.dump(out, open(sys.argv[3], 'w'), indent=1)
