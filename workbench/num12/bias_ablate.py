import numpy as np, sys
from closure import closure
from corrected3 import tree_closure
n, L, s = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
Ws = list(np.load(f"W_n{n}_L{L}_s{s}.npy")); tr = np.load(f"truth_n{n}_L{L}_s{s}.npz"); m = tr["m"][-1]
runs = {"closure": closure(Ws),
        "1pt only": tree_closure(Ws, J=1, rmax=4, two_point=False),
        "2pt only": tree_closure(Ws, J=1, rmax=4, one_point=False),
        "full": tree_closure(Ws, J=1, rmax=4)}
for k, o in runs.items():
    e = o[-1]["m"] - m; mc = o[-1]["m"]
    c = mc @ e/(mc @ mc)
    print(f"{k:9s}: MSE {np.mean(e**2):.2e}  bias {e.mean():+.2e}  scale coef {c:+.4f}  MSE after oracle scale {np.mean((e - c*mc)**2):.2e}")
