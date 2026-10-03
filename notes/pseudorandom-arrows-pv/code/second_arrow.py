# Which second readout kills the dominant second chaos?  A_w = E_X[psi_2(X) h'h'^T]; residual energy (s2^2/2)||A||^2 sin^2(A,B)
import numpy as np
from scipy.stats import norm
rng = np.random.default_rng(2)
n, N = 256, 40000
s2 = 2.0 / n
H = rng.standard_normal((N, n))
for _ in range(4):
    H = np.maximum(H @ (rng.standard_normal((n, n)) * np.sqrt(s2)), 0)
m = H.mean(0); e = m / np.linalg.norm(m); Q = np.eye(n) - np.outer(e, e)
Hp = H @ Q; tau = s2 * (Hp**2).sum(1); st = np.sqrt(tau)
SigP = np.cov(Hp.T, bias=True)                  # Q Sigma Q
M2P = Hp.T @ Hp / N                             # Q M2 Q  (= Q Sigma Q since Qm = 0)
ev, V = np.linalg.eigh(SigP); top = V[:, -1:]
def cos2(A, B): return (np.sum(A * B))**2 / (np.sum(A * A) * np.sum(B * B))
rows = []
for _ in range(30):
    w = rng.standard_normal(n) * np.sqrt(s2)
    a = H @ (np.outer(e, e) @ w)
    psi2 = norm.pdf(a / st) / st
    A = (Hp * psi2[:, None]).T @ Hp / N
    E2 = s2**2 / 2 * np.sum(A * A)
    rows.append([E2, 1 - cos2(A, Q), 1 - cos2(A, SigP)])
r = np.array(rows)
print(f"second-chaos energy (mean over rows)          = {r[:,0].mean():.3e}")
print(f"  remaining after retaining |u|^2 (trace arrow)   : fraction {r[:,1].mean():.3f}")
print(f"  remaining after retaining u^T Sigma' u (variance arrow): fraction {r[:,2].mean():.4f}")
print(f"  Sigma' top eigenvalue {ev[-1]:.2f}, bulk mean {ev[:-1].mean():.3f}")
