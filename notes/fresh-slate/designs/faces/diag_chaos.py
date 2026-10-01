"""How much of each neuron's variance does the second-order chaos in x carry?  |b|^2 (first chaos, face-averaged
gradient norm^2) and 2 tr Q^2 (second chaos: facet births on Gaussian legs) against the closure variance."""
import sys, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "bench")); sys.path.insert(0, HERE)
import bench, fbt
from scipy.special import ndtr
n = int(sys.argv[1]); W = bench.weights_from_seed(11, n, 16).astype(np.float64)
L = 16
mu = np.zeros(n); C = W[0].T @ W[0]; Px = [W[0]]; Pfrom = []; cs = []
for l in range(L):
    if l > 0:
        mu = Ea @ W[l]; C = W[l].T @ Ca @ W[l]
        Px.append((Px[-1] * beta[None, :]) @ W[l]); Pfrom = [(P * beta[None, :]) @ W[l] for P in Pfrom] + [W[l]]; cs.append(cfac)
        # 2 tr Q_j^2 = 2 sum_{(s,k),(s',k')} A_j[sk] A_j[s'k'] (v_sk . v_s'k')^2
        Vs = np.concatenate(Px[:l], 1)                       # (n, l n) facet directions
        A = np.concatenate([Pfrom[s] * cs[s][:, None] for s in range(l)], 0)   # (l n, n) coefficients
        G2 = (Vs.T @ Vs) ** 2
        trQ2 = (A * (G2 @ A)).sum(0)
        b2 = (Px[l] ** 2).sum(0); var = np.diag(C)
        print(f"z-layer {l+1:2d}: first chaos |b|^2/var {np.mean(b2/var):.3f}   second chaos 2trQ^2/var {np.mean(2*trQ2/var):.3f}")
    var = np.diag(C); s_ = np.sqrt(var); t = mu / s_
    Ea = fbt.readout(mu, var, 0 * mu, 0 * mu)
    beta = ndtr(t); cfac = fbt._phi(t) / (2 * s_)
    Ca = fbt.mehler_cov(mu, C, 6)
    Ea2 = (mu * mu + var) * ndtr(t) + mu * s_ * fbt._phi(t); np.fill_diagonal(Ca, Ea2 - Ea * Ea)
