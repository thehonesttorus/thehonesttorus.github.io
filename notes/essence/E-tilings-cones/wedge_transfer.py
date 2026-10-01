"""Parameter-free 'wedge' (two-replica solid angle) model of the transfer coefficients K(l)
measured by the region stream (notes/fresh-slate/breakthrough/region/results/lr_mlp*.json, n = 1024, L = 16).

Conic reading: for a neuron with fresh weight w, P(w.a(x) > 0, w.a(x') > 0) for two independent inputs is the
Gaussian measure of a 2-D wedge of opening pi - arccos(rho) (Sheppard): 1/4 + arcsin(rho)/(2 pi), rho = replica
correlation at that layer (arc-cosine map iterated from rho_0 = 0). Given W, the gates of independent inputs are
independent, so this wedge measure equals the neuron-average of Phi_j^2, and a mean error is transported by
lambda_l = 2 E[Phi^2] = 1/2 + arcsin(rho_l)/pi per layer in energy.  Criticality: rho -> 1, lambda -> 1 - 3/l.
Off-diagonal covariance errors: off->off factor lambda^2, off->diag(next pre-act variance) (8/n)||dC||^2,
diag->mean phi(t)^2/(4 sigma^2) with t ~ N(0, rho/(1-rho)), sigma^2 = 2(1-rho) (He: second moment 2 per neuron).
"""
import json, numpy as np, sys
n, L = 1024, 16
f = lambda r: (np.sqrt(1 - r * r) + (np.pi - np.arccos(r)) * r) / np.pi
rho = [0.0]
for _ in range(L - 1):
    rho.append(f(rho[-1]))
rho = np.array(rho)
lam = 0.5 + np.arcsin(rho) / np.pi
Kmean = np.array([np.prod(lam[s + 1:L]) for s in range(L)])
tau2 = rho / (1 - rho); sig2 = 2 * (1 - rho)
q = (1 / (2 * np.pi)) / np.sqrt(1 + 2 * tau2) / (4 * sig2)   # E[phi(t)^2]/(4 sigma^2)
Koff = np.zeros(L)
for l in range(L - 1):
    tot = 0.0
    for k in range(l + 1, L):
        tot += np.prod(lam[l + 1:k] ** 2) * q[k] * Kmean[k]
    Koff[l] = (8 / n) * tot / n   # final MSE is per neuron: (1/n)||dm||^2
base = '../../fresh-slate/breakthrough/region/results/'
M = [json.load(open(base + f'lr_mlp{i}.json')) for i in range(3)]
print(' l   rho    lambda | Kmean pred  meas(mlp0,1,2)        | Koff pred   meas(mlp0,1,2)')
for l in range(L):
    km = [m['mean'][l] if m['mean'][l] is not None else float('nan') for m in M]
    ko = [m['off'][l] if l < len(m['off']) and m['off'][l] is not None else float('nan') for m in M]
    print(f"{l:2d} {rho[l]:.3f} {lam[l]:.3f} | {Kmean[l]:.3f}  " + ' '.join(f'{v:.3f}' for v in km)
          + f" | {Koff[l]:.2e} " + ' '.join(f'{v:.2e}' for v in ko))
lm = np.log(np.array([[m['mean'][l] for m in M] for l in range(L - 1)]).mean(1)) - np.log(Kmean[:L - 1])
lo = np.log(np.array([[m['off'][l] for m in M] for l in range(L - 1)]).mean(1)) - np.log(Koff[:L - 1])
print('log ratio meas/pred  mean: avg %.2f sd %.2f   off: avg %.2f sd %.2f' % (lm.mean(), lm.std(), lo.mean(), lo.std()))
