"""Profile of the Euler-Stein defect along Eldan's path on a small network (Gaussian closure): g(tau) =
<E_y delta(y), delta(0)> / |delta(0)|^2 for y ~ N(0, tau I), delta(y) = e(y) - y.grad e(y) - Lap e(y) by finite differences.
The weight (1/2)(1+tau)^{-3/2} integrates g to the merge coefficient a; g = (1+tau)^{-1/2} gives a = 1/2.
  python scripts/defect_profile.py n=256 L=4 taus=0.25,1,4,16 samples=3 h=0.02"""
import sys, os, time, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from whest.relu_gauss import relu_moments
kw = dict(n=256, L=4, taus="0.25,1,4,16", samples=3, h=0.02, seed=0)
for a in sys.argv[1:]:
    k, v = a.split("="); kw[k] = type(kw[k])(v)
n, L, h = kw["n"], kw["L"], kw["h"]; taus = [float(x) for x in kw["taus"].split(",")]
rng = np.random.default_rng(kw["seed"]); W = rng.standard_normal((L, n, n)) * np.sqrt(2.0 / n)
def e(y):
    mu = y.copy(); C = np.eye(n)
    for l in range(L):
        mz = W[l] @ mu; Cz = W[l] @ C @ W[l].T; mu, C, _, _ = relu_moments(mz, Cz, 8)
    return mu
def delta(y):
    e0 = e(y); g = np.zeros((n, n)); lap = np.zeros(n)
    for i in range(n):
        d = np.zeros(n); d[i] = h; ep, em = e(y + d), e(y - d); g[i] = (ep - em) / (2 * h); lap += (ep + em - 2 * e0) / h ** 2
    return e0 - y @ g - lap
t0 = time.time(); d0 = delta(np.zeros(n)); print(f"delta(0): rms {np.sqrt(np.mean(d0**2)):.3e} ({time.time()-t0:.0f}s)", flush=True)
print("tau | g(tau) = <E_y delta(y), delta(0)>/|delta(0)|^2 (mean over samples, per-sample values) | |E_y delta|/|delta(0)| | (1+tau)^-1/2")
for tau in taus:
    gs = []; acc = np.zeros(n)
    for s in range(kw["samples"]):
        y = rng.standard_normal(n) * np.sqrt(tau); dy = delta(y); acc += dy; gs.append(np.dot(dy, d0) / np.dot(d0, d0))
    acc /= kw["samples"]
    print(f"{tau:6.2f} | {np.mean(gs):+.3f} ({' '.join(f'{g:+.3f}' for g in gs)}) | {np.sqrt(np.mean(acc**2))/np.sqrt(np.mean(d0**2)):.3f} | {(1+tau)**-0.5:.3f}", flush=True)
