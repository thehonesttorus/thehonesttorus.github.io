# Readout-truncation test: on the SAME samples of the final pre-activation z = W16 h15, compare the sample mean of relu(z)
# with readouts computed from the sample cumulants of z (Gaussian; Edgeworth k3; Edgeworth k3+k4+k3^2).
import numpy as np, time
from scipy.special import ndtr
Wc = np.load("../official/W_off0.npy"); Wf = [np.ascontiguousarray(W.T) for W in Wc]
n = 1024; N = 524288; B = 4096; rng = np.random.default_rng(5)
S = np.zeros((5, n)); R = np.zeros(n); t0 = time.time()
for b in range(N//B):
    h = rng.standard_normal((B, n)).astype(np.float32)
    for W in Wf[:-1]: h = np.maximum(h @ W, 0)
    z = (h @ Wf[-1]).astype(np.float64)
    zp = np.ones_like(z)
    for k in range(1, 5): zp = zp*z; S[k] += zp.sum(0)
    R += np.maximum(z, 0).sum(0)
print(f"forward {time.time()-t0:.0f}s")
m1, m2, m3, m4 = (S[1:]/N); mr = R/N
mu = m1; var = m2 - mu**2; sig = np.sqrt(var)
k3 = m3 - 3*mu*m2 + 2*mu**3; k4 = m4 - 4*mu*m3 + 6*mu**2*m2 - 3*mu**4 - 3*var**2
a = mu/sig; f = np.exp(-a*a/2)/np.sqrt(2*np.pi); P = ndtr(a)
g = sig*(a*P + f)
s3 = k3/(6*sig**3); s4 = k4/(24*sig**4)
e3 = g + s3*(-sig*a*f)
e4 = e3 + s4*(sig*(a*a - 1)*f) + 0.5*s3*s3*(sig*(a**4 - 6*a*a + 3)*f)
for name, r in [("Gaussian (mu, sigma)", g), ("Edgeworth + k3", e3), ("Edgeworth + k3 + k4 + k3^2", e4)]:
    d = r - mr; print(f"{name:30s} readout - sample mean: MSE {np.mean(d**2):.3e}   (scale part {(mr @ d)/(mr @ mr):+.2e})")
print(f"typical standardized cumulants: |k3|/sig^3 median {np.median(np.abs(k3)/sig**3):.3f}, k4/sig^4 median {np.median(k4/sig**4):.4f}, a = mu/sig quantiles {np.quantile(a, [0.05, 0.5, 0.95]).round(2)}")
np.savez("z16_moments_off0.npz", S=S, R=R, N=N)
