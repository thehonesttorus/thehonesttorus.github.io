# Code tree of a finite-memory (de Bruijn) Markov code process, Michon weights kappa(v) = mu(v)^(1/s).
# Claims: zeta(t) = sum_v kappa(v)^t has a simple pole at t = s with residue s/h (h = entropy rate, bits->nats consistent),
# and the subtree residues give the Dixmier measure: Res zeta_v / Res zeta = mu(v).
# Exact: sum_{|v|=k} mu(v)^sigma = pi_sigma^T P_sigma^(k-m) 1 on de Bruijn states (sigma = t/s).
import numpy as np, itertools
rng = np.random.default_rng(4)
m = 3                                   # memory: state = last m bits (de Bruijn graph B(2, m))
S = 2**m
p1 = rng.uniform(0.15, 0.85, S)         # P(next bit = 1 | state)
P = np.zeros((S, S))
for st in range(S):
    for b in (0, 1):
        P[st, ((st << 1) & (S - 1)) | b] = p1[st] if b else 1 - p1[st]
w, V = np.linalg.eig(P.T); pi = np.real(V[:, np.argmax(np.real(w))]); pi /= pi.sum()
h = -sum(pi[st]*(p1[st]*np.log(p1[st]) + (1-p1[st])*np.log(1-p1[st])) for st in range(S))   # nats per symbol
def zeta_sigma(sig, K=20000, start=None):
    # sum over depths k >= m of sum_{|v|=k} mu(v)^sig  (cylinders of length k, first m bits drawn from pi)
    Ps = P**sig; v = (pi if start is None else start)**sig
    tot = 0.0
    for k in range(K):
        tot += v.sum(); v = v @ Ps
    return tot
s = 1.7                                  # free normalisation: spectral dimension
print(f"entropy rate h = {h:.6f} nats/bit;  predicted residue s/h = {s/h:.5f}")
for eps in [0.05, 0.02, 0.01, 0.005]:
    sig = 1 + eps
    z = zeta_sigma(sig)
    print(f"  t = s(1+{eps}): (t - s) * zeta(t) = {s*eps*z:.5f}")
# Dixmier measure of a cylinder v (a de Bruijn state reached at depth m): subtree residue / total residue
eps = 0.002
tot = zeta_sigma(1 + eps)
for st in [0, 3, 5]:
    start = np.zeros(S); start[st] = pi[st]
    sub = zeta_sigma(1 + eps, start=start)
    print(f"  cylinder {st:0{m}b}: Dixmier ratio {sub/tot:.5f}  vs  mu(v) = {pi[st]:.5f}")
