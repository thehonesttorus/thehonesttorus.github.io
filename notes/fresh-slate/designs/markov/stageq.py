"""Stage Q driver: python stageq.py WIDTH SEEDS W_LIST [opts]"""
import sys, time, json, numpy as np
sys.path.insert(0, '.')
from mkv import mkv
from bake import weights
TR = '/tmp/claude-0/-home-user-thehonesttorus-github-io/fc3b2401-9ebc-5116-8952-39e8e61f85e7/scratchpad/truth'
width = int(sys.argv[1]); seeds = [int(s) for s in sys.argv[2].split(',')]
variants = sys.argv[3].split(',')
res = {}
for seed in seeds:
    t = np.load(f'{TR}/w{width}_s{seed}.npz'); W = weights(width, 16, seed)
    noise = t['se2'].mean()
    for var in variants:
        kw = dict(gauss=True, ng_cov=False) if var == 'G' else dict(w=int(var))
        t0 = time.time(); m = mkv(W, **kw); dt = time.time() - t0
        err = ((m - t['means']) ** 2).mean(1)
        res.setdefault(var, []).append(dict(seed=seed, final=err[-1] - noise, layers=err.tolist(), noise=noise, t=dt))
        print(f"w{width} s{seed} {var:>3}: final MSE-noise {err[-1]-noise:.3e} (noise {noise:.1e})  L8 {err[7]:.2e}  L12 {err[11]:.2e}  {dt:.1f}s", flush=True)
for var, r in res.items():
    print(f"== {var}: mean final {np.mean([x['final'] for x in r]):.3e}")
json.dump(res, open(f'results/stageq_w{width}_{"-".join(variants)}.json', 'w'))
