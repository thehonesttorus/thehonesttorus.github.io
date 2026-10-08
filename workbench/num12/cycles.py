# Size of the first cycle (triangle) diagram vs the tree diagrams in kappa_3(z_i), by depth.
import numpy as np, sys
from closure import closure, relu_coeffs
from edgeworth import cumulants_next
from scipy.special import ndtr
n = int(sys.argv[1]); L = int(sys.argv[2]); s = int(sys.argv[3])
Ws = list(np.load(f"W_n{n}_L{L}_s{s}.npy")); o = closure(Ws, keep=True)
rng = np.random.default_rng(0)
for l in [1, 3, 7, 11, 14]:
    if l+1 >= L: continue
    d = o[l]; W = Ws[l+1]; R0 = d["R"].copy(); np.fill_diagonal(R0, 0)
    A = relu_coeffs(d["mu"], d["sig"], 4)
    _, _, parts = cumulants_next(W, d["mu"], d["sig"], d["R"], J=3)
    idx = rng.choice(n, 64, replace=False); tri = []
    for i in idx:
        v = W[i]*A[2]; M = v[:, None]*R0
        tri.append(np.trace(M @ M @ M))
    tri = np.array(tri)
    tree = parts["T3"][idx]; tot = parts["D3"][idx]+parts["P3"][idx]+tree
    print(f"layer {l+1:2d}->{l+2:2d}: rho_max {np.abs(R0).max():.2f} rms rho {np.sqrt(np.mean(R0**2)):.3f} | rms kappa3: all-tree {np.sqrt(np.mean(tot**2)):.2e}, T3 {np.sqrt(np.mean(tree**2)):.2e}, triangle {np.sqrt(np.mean(tri**2)):.2e}  ratio tri/tree-total {np.sqrt(np.mean(tri**2)/np.mean(tot**2)):.2f}")
# depth contraction of the linear-response rows
l = L-1; B = Ws[l]; print("row-norm transport:")
for r in range(min(8, L-1)):
    j = l-1-r; d = o[j]; P = ndtr(d["mu"]/d["sig"])
    nb = (B**2).sum(1)
    pred = 2*((B**2)*(P**2)).sum(1)
    Bn = (B*P) @ Ws[j]
    print(f"  r={r}: mean ||B_i||^2 {nb.mean():.3f} -> next {(Bn**2).sum(1).mean():.3f}  predicted 2<P^2>_B ||B||^2 = {pred.mean():.3f}   2<P(1-P)> = {2*np.mean(P*(1-P)):.3f}")
    B = Bn
