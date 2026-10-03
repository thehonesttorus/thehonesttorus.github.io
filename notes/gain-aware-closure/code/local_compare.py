# Which predicted local defects reproduce the measured local mean defect of each layer?
import numpy as np, sys
from scipy.special import ndtr
from ledger import gstep
from closure import closure
from edgeworth import cumulants_next, edgeworth_shift
n, L, s = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
Ws = list(np.load(f"W_n{n}_L{L}_s{s}.npy")); tr = np.load(f"truth_n{n}_L{L}_s{s}.npz"); hm = np.load(f"hmom_n{n}_L{L}_s{s}.npz")
mstar = tr["m"]; Cstar = [hm["S2"][l] - np.outer(hm["S1"][l], hm["S1"][l]) for l in range(L)]
oc = closure(Ws, keep=True)
print("layer | measured local scale | predicted scale: k3-fresh  k4-fresh  k3+k4 fresh  all(r<=4) | corr(pred_all, measured)")
for l in range(1, L):
    mu, S = Ws[l] @ mstar[l-1], Ws[l] @ Cstar[l-1] @ Ws[l].T
    mloc, _ = gstep(mu, S); dm = mloc - mstar[l]; ms = mstar[l]
    sc = lambda v: (ms @ v)/(ms @ ms)
    sig = np.sqrt(np.diag(S))
    # cumulants of z_l from the closure data of layer l-1 (fresh) and earlier (linear response)
    preds = {}
    k3 = np.zeros(n); k4 = np.zeros(n); B = Ws[l]
    for r in range(min(l, 5)):
        d = oc[l-1-r]; t3, t4, parts = cumulants_next(B, d["mu"], d["sig"], d["R"], J=1)
        if r == 0:
            preds["k3f"] = -edgeworth_shift(mu, sig, t3, 0*t3); preds["k4f"] = -edgeworth_shift(mu, sig, 0*t4, t4)
            preds["f"] = -edgeworth_shift(mu, sig, t3, t4)
        k3 += t3; k4 += t4
        if l-2-r >= 0: B = (B*ndtr(d["mu"]/d["sig"])) @ Ws[l-1-r]
    preds["all"] = -edgeworth_shift(mu, sig, k3, k4)
    c = np.corrcoef(preds["all"], dm)[0, 1]
    print(f" {l+1:3d}  |   {sc(dm):+.5f}           |   {sc(preds['k3f']):+.5f}  {sc(preds['k4f']):+.5f}   {sc(preds['f']):+.5f}    {sc(preds['all']):+.5f} | {c:.2f}")
