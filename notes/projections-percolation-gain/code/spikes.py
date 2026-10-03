# Grok's falsification test: spikes of W^(1) (and of every W) above the Marchenko-Pastur / quarter-circle edge,
# and whether the closure residual lives on outlier directions of the layer covariance.
import numpy as np
from closure import closure
for n, L, s in [(256,16,0),(256,16,1),(256,16,2),(512,16,0),(1024,16,0)]:
    Ws = list(np.load(f"W_n{n}_L{L}_s{s}.npy"))
    edge = 2*np.sqrt(2)                                  # sqrt(n*2/n) * 2 : top singular value edge of N(0,2/n) square matrix
    tops = [np.linalg.svd(W, compute_uv=False)[:3] for W in Ws]
    r = sum(int(np.sum(t > edge*(1 + 3*n**(-2/3)))) for t in tops)
    print(f"n={n} s={s}: edge {edge:.3f}; W1 top singular values {np.round(tops[0],3)}; max over 16 layers {max(t[0] for t in tops):.3f}; spikes above edge*(1+3n^-2/3): {r}")
    tr = np.load(f"truth_n{n}_L{L}_s{s}.npz"); m = tr["m"][-1]
    o = closure(Ws, keep=True)[-1]; mc = o["m"]; C = o["C"]
    e = mc - m; e = e - (mc @ e)/(mc @ mc)*mc       # residual after removing the scale residue
    ev, V = np.linalg.eigh(C); ev, V = ev[::-1], V[:, ::-1]
    for k in [1, 4, 16]:
        frac = np.sum((V[:, :k].T @ e)**2)/np.sum(e**2)
        print(f"     residual energy on top-{k:2d} eigvecs of Cov(h_L): {frac:.3f}  (chance {k/n:.3f}); top eigs/mean eig {np.round(ev[:3]/ev.mean(),1)}")
