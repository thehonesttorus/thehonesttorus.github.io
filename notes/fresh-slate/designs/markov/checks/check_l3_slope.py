"""Layer-3 kappa3 slope vs MC truth for MKV variants (n from argv, truth by MC, 3 layers)."""
import sys, numpy as np
sys.path.insert(0, '.')
import mkv as M
from bake import weights
n = int(sys.argv[1]); N = int(float(sys.argv[2]))
W = weights(n, 16, 1)[:3].astype(np.float64); rng = np.random.default_rng(9)
S = np.zeros((3, n)); cnt = 0
while cnt < N:
    x = rng.standard_normal((100000, n)); x = np.concatenate([x, -x])
    z3 = np.maximum(np.maximum(x @ W[0], 0) @ W[1], 0) @ W[2]
    S[0] += z3.sum(0); S[1] += (z3**2).sum(0); S[2] += (z3**3).sum(0); cnt += len(z3)
m1, m2, m3 = S / cnt; k3 = m3 - 3*m2*m1 + 2*m1**3
for name, kw in [('base', {}), ('ng_site', dict(ng_site=True)), ('no var21', dict(var21=False))]:
    _, dg = M.mkv(W, w=4, var21=kw.pop('var21', 'full'), return_all=True, **kw)
    e = dg[2]['k3']
    print(f"n={n} {name:9s} slope {np.dot(e,k3)/np.dot(k3,k3):.3f}  relerr {np.sqrt(np.mean((e-k3)**2)/np.mean(k3**2)):.3f}", flush=True)
