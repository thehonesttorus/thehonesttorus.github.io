# The "diffraction pattern" of the deep representation: empirical characteristic function of final-layer
# pre-activations along directions, against (i) the moment-matched Gaussian (closure) and (ii) the gain-aware
# closure's scale mixture z = G y, y ~ N(mu, S), G^2 ~ Gamma(1/gamma, gamma), with gamma from the weights-only formula.
import numpy as np
from math import lgamma
from gac import gac
n, L, s, T = 256, 16, 0, 400000
Ws = np.load(f"W_n{n}_L{L}_s{s}.npy"); Wf = [W.T.astype(np.float32) for W in Ws]
g = gac(list(Ws))[-1]; mu, S, gam = g["mu"], g["S"], g["gamma"]
print(f"GAC gain variance at the final pre-activation: gamma = {gam:.4f}")
rng = np.random.default_rng(5); Z = []
for _ in range(T//20000):
    h = rng.standard_normal((20000, n)).astype(np.float32)
    for W in Wf[:-1]: h = np.maximum(h @ W, 0)
    Z.append((h @ Wf[-1]).astype(np.float64))
Z = np.concatenate(Z)
# Gamma(k, 1/k) quadrature for G^2 (generalised Laguerre)
k = 1/gam; xg, wg = np.polynomial.laguerre.laggauss(80)
# E f(G^2) = int x^{k-1} e^{-kx} k^k/Gamma(k) f(x) dx ; substitute y = kx
w = wg*np.exp((k-1)*np.log(xg/k) - lgamma(k) + np.log(k) - np.log(k)) ; w = w/ w.sum(); G2 = xg/k
ts = np.array([1.0, 2.0, 2.5, 3.0, 3.5, 4.0])
rows = {"empirical": [], "Gaussian": [], "GAC mixture": []}
dirs = [rng.standard_normal(n) for _ in range(60)] + [np.eye(n)[i] for i in range(0, n, 4)]
for v in dirs:
    y = Z @ v; m, sd = y.mean(), y.std()
    rows["empirical"].append([abs(np.mean(np.exp(1j*t*(y-m)/sd))) for t in ts])
    rows["Gaussian"].append(np.exp(-ts**2/2))
    a, b2 = v @ mu, v @ S @ v; G = np.sqrt(G2)
    rows["GAC mixture"].append([abs(np.sum(w*np.exp(1j*t*G*a/sd - t*t*G2*b2/(2*sd*sd))*np.exp(-1j*t*m/sd))) for t in ts])
print(f"T={T}, noise floor ~{1/np.sqrt(T):.4f}; mean over {len(dirs)} directions (60 random, {len(dirs)-60} single neurons)")
print("      t     " + "  ".join(f"{t:7.1f}" for t in ts))
for k_, r in rows.items(): print(f"{k_:12s}" + "  ".join(f"{x:7.4f}" for x in np.mean(r, 0)))
e = np.array(rows["empirical"]); 
for k_ in ["Gaussian", "GAC mixture"]:
    print(f"rms log-error at t=2..3.5 vs empirical: {k_:12s} {np.sqrt(np.mean(np.log(np.array(rows[k_])[:, 1:5]/e[:, 1:5])**2)):.3f}")
