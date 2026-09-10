# Requirement 4, with the actual submitted estimator (numpy mode, identical equations to the
# flopscope path) over all 100 mini-split networks against their 1e9-sample ground truth.
import numpy as np, time, sys
sys.path.insert(0,'.')
import estimator as E
from types import SimpleNamespace
B=np.load('../bank.npz'); seeds=B['seeds']; Y=B['Y'].astype(np.float64); names=B['names']
n,L=1024,16
def regen(seed):
    ss=np.random.SeedSequence(int(seed)).spawn(3); rng=np.random.default_rng(ss[0]); sc=float(np.sqrt(2.0/n))
    return [(rng.standard_normal((n,n))*sc).astype(np.float32) for _ in range(L)]
def score(kernel):
    v=[]; t0=time.time()
    for k in range(len(seeds)):
        W=regen(seeds[k])
        mlp=SimpleNamespace(width=n,depth=L,weights=W,seed=int(seeds[k])&0xFFFFFFFF,name=str(names[k]))
        est=E.FiniteResolutionCumulantEstimator(xp=np); est.USE_KERNEL=kernel
        P=est.predict(mlp,2**41)
        v.append(np.mean((P.astype(np.float64)-Y[k])**2,axis=1)[-1])
        if k%25==24: print(f"    {k+1}/100  {time.time()-t0:.0f}s", flush=True)
    return np.array(v)
for kernel in (False, True):
    v=score(kernel)
    tag="with shipped kernel" if kernel else "without kernel (previous)"
    print(f"{tag:26s} mean {v.mean():.4e}  median {np.median(v):.4e}  sd {v.std():.3e}"
          f"  fitted(0-7) {v[:8].mean():.4e}  held-out(8-99) {v[8:].mean():.4e}", flush=True)
    np.save(f'bank_est_{int(kernel)}.npy', v)
