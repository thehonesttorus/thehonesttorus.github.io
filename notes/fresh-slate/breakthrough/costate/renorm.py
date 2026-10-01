"""Is the old content's readout effect a rescaling of something cheap? Leave-one-out ensemble-fitted scalar gamma."""
import sys, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "../../bench"))
import bench
P = os.environ["PRED_DIR"]
name = sys.argv[1]
S = bench.load_set(name); M = len(S["seeds"])
ld = lambda v, i: np.load(os.path.join(P, f"{name}_{v}_{i}.npy"))[-1]
mlps = [i for i in range(M) if os.path.exists(os.path.join(P, f"{name}_full_{i}.npy"))]
print(f"{name}: {len(mlps)} MLPs")
for A in ("A0", "A1", "A3"):
    for cheap in (A + "slice", "young"):
        rows = []
        for i in mlps:
            t = S["means"][i][-1]; full = ld("full", i); pa = ld(A, i); g = ld("gauss", i)
            O = full - pa
            X = (ld(cheap, i) - pa) if cheap != "young" else (pa - g)
            rows.append((t, full, pa, O, X))
        cs = [np.dot(O, X) / np.linalg.norm(O) / np.linalg.norm(X) for _, _, _, O, X in rows]
        mse_A = np.mean([np.mean((pa - t) ** 2) for t, _, pa, _, _ in rows])
        mse_full = np.mean([np.mean((f - t) ** 2) for t, f, _, _, _ in rows])
        loo = []
        for k in range(len(rows)):
            num = sum(np.dot(r[3], r[4]) for j, r in enumerate(rows) if j != k)
            den = sum(np.dot(r[4], r[4]) for j, r in enumerate(rows) if j != k)
            gam = num / den
            t, f, pa, O, X = rows[k]
            loo.append(np.mean((pa + gam * X - t) ** 2))
        print(f"  {A} + gamma*({cheap:8s}): cos(old, cheap) = {np.mean(cs):+.3f}   mse {A} {mse_A:.3e} -> LOO-renormalised {np.mean(loo):.3e}  (all pairs {mse_full:.3e})")
