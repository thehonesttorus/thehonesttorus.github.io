# Is the measured local mean defect the scale-mixture term (gamma_l/8) sigma phi(a)(1+a^2) with the ACCUMULATED gain
# variance gamma_l (measured), plus the tree-level (fresh + linear response) part?
import numpy as np
from scipy.special import ndtr
from ledger import gstep
from closure import closure, phi
from edgeworth import cumulants_next, edgeworth_shift
n, L, s = 256, 16, 0
gam = [0, 0.0156, 0.0258, 0.0331, 0.0395, 0.0432, 0.0492, 0.0462, 0.0485, 0.0520, 0.0526, 0.0519, 0.0573, 0.0551, 0.0577, 0.0552]
Ws = list(np.load(f"W_n{n}_L{L}_s{s}.npy")); tr = np.load(f"truth_n{n}_L{L}_s{s}.npz"); hm = np.load(f"hmom_n{n}_L{L}_s{s}.npz")
mstar = tr["m"]; Cstar = [hm["S2"][l] - np.outer(hm["S1"][l], hm["S1"][l]) for l in range(L)]
oc = closure(Ws, keep=True)
print("layer | measured | mixture(gamma_l)  tree(r<=4)  mixture+tree(prev-gamma) | corr: mixture  tree  sum")
for l in range(1, L):
    mu, S = Ws[l] @ mstar[l-1], Ws[l] @ Cstar[l-1] @ Ws[l].T
    mloc, _ = gstep(mu, S); dm = mloc - mstar[l]; ms = mstar[l]; sc = lambda v: (ms @ v)/(ms @ ms)
    sig = np.sqrt(np.diag(S)); a = mu/sig
    mix = gam[l]/8*sig*phi(a)*(1+a*a)
    mixp = gam[l-1]/8*sig*phi(a)*(1+a*a)        # gain accumulated up to the previous layer (fresh part from trees)
    k3 = np.zeros(n); k4 = np.zeros(n); B = Ws[l]
    for r in range(min(l, 5)):
        d = oc[l-1-r]; t3, t4, _ = cumulants_next(B, d["mu"], d["sig"], d["R"], J=1); k3 += t3; k4 += t4
        if l-2-r >= 0: B = (B*ndtr(d["mu"]/d["sig"])) @ Ws[l-1-r]
    tree = -edgeworth_shift(mu, sig, k3, k4)
    cc = lambda p: np.corrcoef(p, dm)[0, 1]
    print(f" {l+1:3d}  | {sc(dm):+.5f} |  {sc(mix):+.5f}        {sc(tree):+.5f}     {sc(mixp+tree):+.5f}            | {cc(mix):.2f}  {cc(tree):.2f}  {cc(mixp+tree):.2f}")
