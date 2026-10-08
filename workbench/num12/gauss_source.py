# Is the tree-level cumulant formula accurate when the source layer is exactly Gaussian but strongly correlated?
import numpy as np, sys
from closure import closure, relu_var
from edgeworth import cumulants_next, edgeworth_shift
from twopoint import twopoint
n = 256; Ws = list(np.load("W_n256_L16_s0.npy")); o = closure(Ws, keep=True)
for l in [1, 5, 11]:
    d = o[l]; W = Ws[l+1]; L = np.linalg.cholesky(d["mu"]*0 + (d["R"]*np.outer(d["sig"], d["sig"])) + 1e-12*np.eye(n))
    rng = np.random.default_rng(1); T = 4_000_000; B = 20000
    S = np.zeros((4, n)); G = np.zeros(n)
    mu_z = W @ d["m"]
    for _ in range(T//B):
        y = d["mu"] + rng.standard_normal((B, n)) @ L.T
        z = np.maximum(y, 0) @ W.T - mu_z
        for p in range(4): S[p] += (z**(p+1)).sum(0)
        G += np.maximum(z + mu_z, 0).sum(0)
    S /= T; G /= T
    var = S[1]-S[0]**2; k3 = S[2]-3*S[0]*S[1]+2*S[0]**3
    k4 = S[3]-4*S[0]*S[2]+6*S[0]**2*S[1]-3*S[0]**4-3*var**2
    p3, p4, parts = cumulants_next(W, d["mu"], d["sig"], d["R"], J=4)
    sig_z = np.sqrt(np.diag(W @ d["C"] @ W.T))
    from closure import relu_coeffs
    M0 = relu_coeffs(mu_z, sig_z, 2)[0]
    dm_pred = edgeworth_shift(mu_z, sig_z, p3, p4)
    print(f"source layer {l+1} (R rms {np.sqrt(np.mean((d['R']-np.eye(n))**2)):.3f}):")
    print(f"   k3: true rms {np.sqrt(np.mean(k3**2)):.3e}  pred err {np.sqrt(np.mean((p3-k3)**2)):.3e}  noise ~{np.sqrt(np.mean(6*var**3/T)):.1e}")
    print(f"   k4: true mean {k4.mean():.3e}  pred mean {p4.mean():.3e}  err rms {np.sqrt(np.mean((p4-k4)**2)):.3e}  noise ~{np.sqrt(np.mean(96*var**4/T)):.1e}")
    print(f"   mean shift: true rms {np.sqrt(np.mean((G-M0)**2)):.3e}  pred err {np.sqrt(np.mean((G-M0-dm_pred)**2)):.3e}  noise ~{np.sqrt(np.mean(var/T)):.1e}")
