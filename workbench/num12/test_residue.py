import numpy as np, sys
from residue import closure_with_residue
n, L, s = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
Ws = list(np.load(f"W_n{n}_L{L}_s{s}.npy")); tr = np.load(f"truth_n{n}_L{L}_s{s}.npz")
o = closure_with_residue(Ws)
print(f"n={n} L={L} seed={s}")
for l in range(L):
    m = tr["m"][l]; mc = o[l]["m_closure"]; c = mc @ (mc-m)/(mc @ mc)
    if l in (1, 3, 7, 11, L-1) or l == L-1:
        print(f"  layer {l+1:2d}: oracle c {c:+.4f}  gamma/8 from weights {o[l]['gamma']/8:+.4f} | MSE closure {np.mean((mc-m)**2):.2e} -> residue-corrected {np.mean((o[l]['m']-m)**2):.2e} (oracle scale {np.mean((mc*(1-c)-m)**2):.2e}, MC@B {tr['v'][l].mean()/65536:.2e})")
