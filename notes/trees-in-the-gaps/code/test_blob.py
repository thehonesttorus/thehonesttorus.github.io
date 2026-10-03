import numpy as np
from scipy.special import ndtr
from closure import closure, phi
from edgeworth import cumulants_next
from blob import g1_kappa3
n, L, s = 256, 4, 0
Ws = list(np.load(f"W_n{n}_L{L}_s{s}.npy")); z = np.load(f"zstats_n{n}_L{L}_s{s}.npz")
oc = closure(Ws, keep=True)
for l in range(2, L):
    k3 = np.zeros(n); B = Ws[l]
    for r in range(l):
        d = oc[l-1-r]; t3, t4, _ = cumulants_next(B, d["mu"], d["sig"], d["R"]); k3 += t3
        if l-2-r >= 0: B = (B*ndtr(d["mu"]/d["sig"])) @ Ws[l-1-r]
    d1 = oc[l-1]; d0 = oc[l-2]
    P = ndtr(d1["mu"]/d1["sig"]); p = phi(d1["mu"]/d1["sig"])/d1["sig"]
    g1 = g1_kappa3(Ws[l], Ws[l-1], P, p, d0["mu"], d0["sig"], d0["R"])
    tr = z["k3"][l]
    print(f"layer {l+1}: true rms {np.sqrt(np.mean(tr**2)):.3e} | err trees+LR {np.sqrt(np.mean((k3-tr)**2)):.3e} | G1 rms {np.sqrt(np.mean(g1**2)):.3e} | err +G1 {np.sqrt(np.mean((k3+g1-tr)**2)):.3e} | corr(G1, residual) {np.corrcoef(g1, tr-k3)[0,1]:.2f}")
