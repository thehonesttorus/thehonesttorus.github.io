# The chain's regenerated kappa4 off-diagonal core G_off = LAM[l] C_off (LAM fitted offline): what is its physical value?
# Fit the TRUE pair fourth cumulant of h_l, C4_ac = kappa(h_a,h_a,h_c,h_c), against Cov(h_a,h_c) off the diagonal, and
# decompose C4 into the Gaussian one-loop part, the kappa3(y)-variation (scale-mixture slices at the measured g_l, and the
# tree remainder) and the kappa4(y)-variation.
import numpy as np, sys, time
sys.path.insert(0, "../num12")
from closure import relu_coeffs, relu2_coeffs
from pairvar import Layer, y_cumulants, h_cumulants, sm_slices
net = int(sys.argv[1]) if len(sys.argv) > 1 else 0
D = np.load(f"mc_cum_off{net}.npz"); L, n = D["s1y"].shape
gm = np.load(f"k3channel_off{net}.npz")["gm"]
LAM = [1.9516e-03, 7.5412e-03, 9.9784e-03, 1.0922e-02, 1.1267e-02, 1.1241e-02, 1.0951e-02, 1.0348e-02, 9.8982e-03, 9.3751e-03, 8.9429e-03, 8.5313e-03, 8.1677e-03, 7.7712e-03, 7.3375e-03, 6.9496e-03]
iu = np.triu_indices(n, 1)
def fit(x, y): return (x @ y) / (x @ x)
print(f"net {net}: pair kappa4 of h_l against Cov(h_a,h_c): lambda = <C4, C>/<C, C> off-diagonal, and R^2")
print(" l | lambda true (R^2) | one-loop Gaussian | kappa3 var: scale-mixture(g_l) / tree | kappa4 var | sum | chain LAM[l] | g_l")
t0 = time.time()
for l in range(L):
    Y = y_cumulants(D, l); mu, S = Y["mu"], Y["S"]; sig = np.sqrt(np.diag(S))
    A = relu_coeffs(mu, sig, 19); B = relu2_coeffs(mu, sig, 19); lay = Layer(mu, S, A, B)
    C3t, C4t, k4t = h_cumulants(D, l); Ch = D["S2h"][l] - np.outer(D["s1h"][l], D["s1h"][l]); cv = Ch[iu]
    lt = fit(cv, C4t[iu]); r2 = 1 - np.var(C4t[iu] - lt * cv) / np.var(C4t[iu])
    l1 = fit(cv, lay.C4G[iu])
    sl = sm_slices(mu, S, gm[l]); _, d4g, _ = lay.var3(sl["k3"], sl["K21"]); _, d4t, _ = lay.var3(Y["k3"] - sl["k3"], Y["K21"] - sl["K21"]); _, d44, _ = lay.var4(Y["k4"], Y["K31"], Y["K22"])
    lg, ltr, l4 = fit(cv, d4g[iu]), fit(cv, d4t[iu]), fit(cv, d44[iu])
    print(f" {l:2d} | {lt:.5f} ({r2:.3f}) | {l1:+.5f} | {lg:+.5f} / {ltr:+.5f} | {l4:+.5f} | {l1+lg+ltr+l4:.5f} | {LAM[l]:.5f} | {gm[l]:.5f}   [{time.time()-t0:.0f}s]", flush=True)
