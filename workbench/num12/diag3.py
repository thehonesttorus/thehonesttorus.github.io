import numpy as np
from math import factorial
from scipy.special import ndtr
from closure import relu_coeffs, relu_var
from edgeworth import cumulants_next
from corrected2 import corrected_closure2
n, L, s = 256, 4, 0
Ws = list(np.load(f"W_n{n}_L{L}_s{s}.npy")); z = np.load(f"zstats_n{n}_L{L}_s{s}.npz")
o = corrected_closure2(Ws)
for l in range(L):
    dmu = o[l]["mu"]-z["mu"][l]; dsig2 = o[l]["sig"]**2-z["var"][l]
    print(f"layer {l+1}: mu err rms {np.sqrt(np.mean(dmu**2)):.2e} | var err rms {np.sqrt(np.mean(dsig2**2)):.2e} mean {dsig2.mean():+.2e} (var~{z['var'][l].mean():.3f})")
# predicted cumulants with response at each layer
hist = []
from closure import closure
oc = closure(Ws, keep=True)
for l in range(1, L):
    k3 = np.zeros(n); k4 = np.zeros(n); B = Ws[l]
    for r in range(l):
        d = oc[l-1-r]; t3, t4, _ = cumulants_next(B, d["mu"], d["sig"], d["R"]); k3 += t3; k4 += t4
        if l-2-r >= 0: B = (B*ndtr(d["mu"]/d["sig"])) @ Ws[l-1-r]
        if r == 0: k30, k40 = k3.copy(), k4.copy()
    T = z["T"]; vz = z["var"][l]
    n3 = np.sqrt(np.mean(6*vz**3/T)); n4 = np.sqrt(np.mean(96*vz**4/T))
    print(f"layer {l+1}: k3 true rms {np.sqrt(np.mean(z['k3'][l]**2)):.2e} err(r=0) {np.sqrt(np.mean((k30-z['k3'][l])**2)):.2e} err(all r) {np.sqrt(np.mean((k3-z['k3'][l])**2)):.2e} [noise {n3:.1e}]")
    print(f"          k4 true mean {z['k4'][l].mean():.2e} err(r=0) {np.sqrt(np.mean((k40-z['k4'][l])**2)):.2e} err(all r) {np.sqrt(np.mean((k4-z['k4'][l])**2)):.2e} [noise {n4:.1e}]")
