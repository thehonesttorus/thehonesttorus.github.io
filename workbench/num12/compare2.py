import numpy as np, sys
from closure import closure
from corrected import corrected_closure
n, L, s = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
Ws = list(np.load(f"W_n{n}_L{L}_s{s}.npy")); tr = np.load(f"truth_n{n}_L{L}_s{s}.npz")
m, v, T = tr["m"], tr["v"], int(tr["T"])
o0 = closure(Ws); o1 = corrected_closure(Ws); o2 = corrected_closure(Ws, var_corr=False)
for l in range(L):
    mc = v[l].mean()/65536
    f = lambda o: np.mean((o[l]["m"]-m[l])**2)
    print(f"layer {l+1:2d}: MC@B {mc:.2e} | closure {f(o0):.2e} | +tree (mean only) {f(o2):.2e} | +tree (mean+var) {f(o1):.2e} | noise {v[l].mean()/T:.1e}")
