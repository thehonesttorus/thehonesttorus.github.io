"""Stage Q over own truths: python stageq2.py WIDTH SEEDS VARIANTS -> results/q_w{W}.json"""
import sys, time, json, os, numpy as np
sys.path.insert(0, '.')
from mkv import mkv
from mkv2 import mkv2
from bake import weights
TR = '/tmp/claude-0/-home-user-thehonesttorus-github-io/fc3b2401-9ebc-5116-8952-39e8e61f85e7/scratchpad/truth'
VAR = {'G': dict(gauss=True, ng_cov=False, var21=False), 'w1': dict(w=1, var21='full'),
       'w2': dict(w=2, var21='full'), 'w4': dict(w=4, var21='full'), 'w16': dict(w=16, var21='full'),
       'w16n2': dict(w=16, var21='full', two_site=False),
       'c1': dict(w=1, v2=1), 'c2': dict(w=2, v2=1), 'c4': dict(w=4, v2=1), 'c16': dict(w=16, v2=1)}
width = int(sys.argv[1]); seeds = [int(s) for s in sys.argv[2].split(',')]; vs = sys.argv[3].split(',')
out = f'results/q_w{width}.json'
res = json.load(open(out)) if os.path.exists(out) else {}
for seed in seeds:
    t = np.load(f'{TR}/w{width}_s{seed}.npz'); W = weights(width, 16, seed); nz = float(t['se2'].mean())
    for v in vs:
        t0 = time.time(); kw = dict(VAR[v]); m = mkv2(W, w=kw['w']) if kw.pop('v2', 0) else mkv(W, **kw); e = ((m - t['means']) ** 2).mean(1)
        res.setdefault(v, {})[str(seed)] = dict(raw=float(e[-1] - nz), noise=nz, layers=e.tolist(), t=time.time() - t0)
        print(f"w{width} s{seed} {v:6s} raw {e[-1]-nz:.3e} noise {nz:.1e} all-layer {e.mean():.2e} {time.time()-t0:.0f}s", flush=True)
        json.dump(res, open(out, 'w'), indent=1)
