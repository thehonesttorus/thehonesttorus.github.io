import numpy as np, time
B=np.load('../bank.npz'); seeds=B['seeds']; Y=B['Y']
import twin2 as T
T.load_nets=lambda: (list(seeds), Y.astype(np.float64))
import ablate as AB, importlib; importlib.reload(AB)
import cov_rho as CR; importlib.reload(CR)
t0=time.time()
for name,kw in (("A4 only (run19)",dict(S3=0.0)), ("A4 + S3=0.22 (run20)",dict(S3=0.22))):
    v=np.array([CR.predict(k,**kw) for k in range(100)])
    print(f"{name:24s} mean {v.mean():.4e}  median {np.median(v):.4e}"
          f"  fitted(0,2,4,6) {v[[0,2,4,6]].mean():.4e}  held out(8-99) {v[8:].mean():.4e}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
