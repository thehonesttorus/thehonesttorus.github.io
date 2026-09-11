# Full 100-network mini split, 1e9-sample ground truth: the shipped rank-one fourth-cumulant term.
import numpy as np, time
B=np.load('../bank.npz'); seeds=B['seeds']; Y=B['Y']
import twin2 as T
T.load_nets=lambda: (list(seeds), Y.astype(np.float64))
import ablate as AB, importlib; importlib.reload(AB)
import cov_terms as CT; importlib.reload(CT)
t0=time.time()
for name,co in (("deployed (no A4 term)",{}), ("with A4 = 0.0020",{'T1':0.0020})):
    v=np.array([CT.predict(k,co) for k in range(100)])
    print(f"{name:24s} mean {v.mean():.4e}  median {np.median(v):.4e}"
          f"  fitted(0,2,4,6) {v[[0,2,4,6]].mean():.4e}  held-out(8-99) {v[8:].mean():.4e}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
