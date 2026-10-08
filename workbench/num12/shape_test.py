import numpy as np
from residue import closure_with_residue
from closure import phi
cases = [(256,4,0),(256,8,0),(256,16,0),(256,16,1),(512,16,0)]
print("case            closure   uniform(t=1)  gap(t=1)   uniform(.95)  gap(.95)   oracle-uniform  oracle-gap")
for n, L, s in cases:
    Ws = list(np.load(f"W_n{n}_L{L}_s{s}.npy")); tr = np.load(f"truth_n{n}_L{L}_s{s}.npz"); m = tr["m"][-1]
    row = []
    for tau in [1.0, 0.95]:
        o = closure_with_residue(Ws, tau=tau)[-1]; g = o["gamma"]; mc = o["m_closure"]; a = o["mu"]/o["sig"]
        shape = o["sig"]*phi(a)*(1+a*a)
        row += [np.mean((mc*(1-g/8)-m)**2), np.mean((mc - g/8*shape - m)**2)]
    e = mc - m; cu = mc @ e/(mc @ mc); cg = shape @ e/(shape @ shape)
    print(f"n={n} L={L:2d} s={s}  {np.mean(e**2):.2e}  " + "  ".join(f"{x:.2e}  " for x in row) + f"  {np.mean((e-cu*mc)**2):.2e}     {np.mean((e-cg*shape)**2):.2e}")
