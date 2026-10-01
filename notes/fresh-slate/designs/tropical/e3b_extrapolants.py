"""E3b: how much of the T = 0 closure error can extrapolation of closure values G(T) (T > 0) remove?
Only G (cheap) is needed; the target is the bench truth at T = 0.  Families: polynomial in T^2 or T on a window
[T0, T1], several degrees; diagonal Pade in T.  Reports mean final-layer MSE over a set's MLPs per extrapolant, and the
best one (selection on the same truth: an optimistic bound)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../bench"))
import numpy as np, bench
from e3_temperature import closure
Ts = np.round(np.concatenate([[0.0], np.arange(0.025, 1.21, 0.025)]), 4)
res = {}
for setname in sys.argv[1].split(","):
    S = bench.load_set(setname)
    for i in range(len(S["seeds"])):
        Ws = [w.astype(np.float64) for w in bench.weights(S, i)]
        G = np.array([closure(Ws, T)[-1] for T in Ts]); t = S["means"][i][-1]
        mse = lambda p: ((p - t) ** 2).mean()
        res.setdefault((setname, "closure T=0"), []).append(mse(G[0]))
        for var in ("T2", "T"):
            x = Ts ** 2 if var == "T2" else Ts
            for T0 in (0.025, 0.1, 0.2, 0.3):
                for T1 in (0.4, 0.6, 0.9, 1.2):
                    idx = (Ts >= T0) & (Ts <= T1)
                    for deg in (1, 2, 3, 4, 6):
                        if idx.sum() <= deg + 2: continue
                        c = np.polynomial.polynomial.polyfit(x[idx] - x[idx].mean(), G[idx], deg)
                        p = np.polynomial.polynomial.polyval(-x[idx].mean(), c)
                        res.setdefault((setname, f"poly {var} [{T0},{T1}] d{deg}"), []).append(mse(p))
    print(setname, "closure T=0: %.3e" % np.mean(res[(setname, "closure T=0")]))
    best = sorted([(np.mean(v), k[1]) for k, v in res.items() if k[0] == setname])[:6]
    for v, k in best: print("   %.3e  %s" % (v, k))
