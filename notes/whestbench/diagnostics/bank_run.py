# Requirement 4: held-out bank of all 100 mini-split networks, 1e9-sample ground truth
# (floor 7.5e-11, so 1e-7-scale differences are resolved ~1000x over).
import numpy as np, time, sys
import ladder
B=np.load('../bank.npz'); seeds=B['seeds']; Y=B['Y']; names=B['names']
import twin2 as T
T.load_nets = lambda: (list(seeds), Y.astype(np.float64))
import importlib; importlib.reload(ladder)
cfgs=[("official covariance propagation (reference)", None),
      ("plain K2 Gaussian closure",              dict(K=0,kernel=False,mc=False)),
      ("K2 + shipped rank-1 kernel",             dict(K=0,kernel=True, mc=False)),
      ("age-1 source + shipped kernel",          dict(K=1,kernel=True, mc=False)),
      ("deployed: age-1 + kernel + MC channels", dict(K=1,kernel=True, mc=True,N_MC=2400))]
res={}
t0=time.time()
for name,kw in cfgs:
    if kw is None: continue
    v=np.array([ladder.run(k,**kw) for k in range(100)])
    res[name]=v
    print(f"{name:42s} mean {v.mean():.3e}  median {np.median(v):.3e}  sd {v.std():.3e}"
          f"  fitted(0-7) {v[:8].mean():.3e}  held-out(8-99) {v[8:].mean():.3e}  [{time.time()-t0:.0f}s]", flush=True)
np.savez('bank_results.npz', **{k:v for k,v in res.items()})
