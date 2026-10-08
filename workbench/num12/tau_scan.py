import numpy as np
from residue import closure_with_residue
cases = [(256,4,0),(256,8,0),(256,16,0),(256,16,1),(512,16,0)]
print("case          closure    tau=1.00   tau=0.95   tau=0.90   oracle-scale  MC@B")
for n, L, s in cases:
    Ws = list(np.load(f"W_n{n}_L{L}_s{s}.npy")); tr = np.load(f"truth_n{n}_L{L}_s{s}.npz"); m = tr["m"][-1]
    row = []
    for tau in [1.0, 0.95, 0.90]:
        o = closure_with_residue(Ws, tau=tau)[-1]; row.append(np.mean((o["m"]-m)**2))
    mc = o["m_closure"]; c = mc @ (mc-m)/(mc @ mc)
    print(f"n={n} L={L:2d} s={s}  {np.mean((mc-m)**2):.2e}   " + "   ".join(f"{x:.2e}" for x in row) + f"   {np.mean((mc*(1-c)-m)**2):.2e}   {tr['v'][-1].mean()/65536:.2e}")
