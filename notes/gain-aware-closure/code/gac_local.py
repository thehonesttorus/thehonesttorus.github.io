# Local defect of one step, same footing as bundle_local.py: closure step vs gain-aware step (1-D base = the gain),
# from the true moments of h_l (hmoments, T=4e6) against the true next-layer mean (truth, T=1.6e7).
import numpy as np
from scipy.special import ndtr
from closure import phi
from gac import gac, EG
n, L, s = 256, 16, 0
Ws = np.load(f"W_n{n}_L{L}_s{s}.npy"); hm = np.load(f"hmom_n{n}_L{L}_s{s}.npz"); mt = np.load(f"truth_n{n}_L{L}_s{s}.npz")["m"]
gmeas = [0, .0156, .0258, .0331, .0395, .0432, .0492, .0462, .0485, .0520, .0526, .0519, .0573, .0551, .0577, .0552]
gform = [d["gamma"] for d in gac(list(Ws))]
def mstep(mu, var): sig = np.sqrt(var); a = mu/sig; return sig*(a*ndtr(a) + phi(a))
print(" layer | closure rms 8*scale | GAC-step (measured gamma): MSE/closure 8*scale | GAC-step (weights-only gamma): MSE/closure 8*scale")
for l in [0, 2, 4, 6, 10, 14]:
    W = Ws[l+1]; m = hm["S1"][l]; C = hm["S2"][l] - np.outer(m, m); ms = mt[l+1]; sc = lambda d: 8*(ms @ d)/(ms @ ms)
    nu = W @ m; v = np.einsum("ij,jk,ik->i", W, C, W)
    d0 = mstep(nu, v) - ms; row = f"  {l+2:3d}  | {np.sqrt(np.mean(d0**2)):.2e} {sc(d0):+.4f} |"
    for g in [gmeas[l+1], gform[l+1]]:
        g1 = EG(g); mu = nu/g1; d = g1*mstep(mu, v + nu**2 - mu**2) - ms
        row += f"        {np.mean(d**2)/np.mean(d0**2):.3f}    {sc(d):+.4f}          |"
    print(row)
