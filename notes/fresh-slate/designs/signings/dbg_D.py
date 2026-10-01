import numpy as np, stageq, copula1 as c1
from copula import FACT, PQ
W, mean, noise, _ = stageq.load(64, 0); n = 64
d = np.load("/root/sg/mcstats_w64s0_4000000.npz"); S, C, C21, done = d['S'], d['C'], d['C21'], int(d['done'])
l = 1; m = S[l]/done; mu = m[0]; E2 = C[l]/done; E21 = C21[l]/done
Dmc = E21 - np.outer(m[1], mu) - 2*E2*mu[None,:] + 2*np.outer(mu**2, mu)
off = ~np.eye(n, dtype=bool); r = lambda x: np.sqrt(np.mean(x[off]**2))
C1 = W[0].T@W[0]; s = np.sqrt(np.diag(C1)); R = C1/np.outer(s,s)
pr = c1.ext_profiles(np.zeros(n), s, np.tile([1.,0,0],(n,1)))
h, F2, F3 = pr['h'], pr['F2'], pr['F3']; mu_a = h[:,0]
A2 = F2 - 2*mu_a[:,None]*h; A2[:,0] += mu_a**2
A3 = F3 - 3*mu_a[:,None]*F2 + 3*(mu_a**2)[:,None]*h; A3[:,0] -= mu_a**3
k3a = A3[:,0]; R0 = R.copy(); np.fill_diagonal(R0,0)
w = W[1]; W2 = w**2
K21 = c1.M(A2, h, R0)*off; K21W = K21@w
parts = {}
parts['aaa'] = W2.T@(k3a[:,None]*w)
parts['aab'] = W2.T@K21W
parts['aca'] = 2*(w*K21W).T@w
X = {p: (R0**p)@(w*h[:,p:p+1])/FACT[p] for p in range(1,PQ+1)}
pc = 0; pa = 0
for p in range(1,PQ+1):
    for q in range(1,PQ+1):
        pc = pc + (h[:,p+q][:,None]*X[p]*X[q]).T@w
        pa = pa + 2*(w*h[:,p+q][:,None]*X[p]).T@X[q]
parts['pathC'] = pc; parts['pathA'] = pa
tot = 0
print("MC rms", r(Dmc))
for k,v in parts.items():
    tot = tot + v; print(f"{k:6s} rms {r(v):.3e}  cum err {r(tot-Dmc):.3e}")
zc = 0; za = 0
for p in range(1,PQ+1):
    for q in range(1,PQ+1):
        Z = (R0**(p+q))@(W2*(h[:,p]*h[:,q])[:,None])/(FACT[p]*FACT[q])
        zc = zc + (h[:,p+q][:,None]*Z).T@w
        Xpq = (R0**(p+q))@(w*h[:,p+q][:,None])            # sum_a rho_ab^{p+q} W_aj h_a(p+q), rows b
        za = za + 2*(w*(h[:,p]*h[:,q])[:,None]*Xpq/(FACT[p]*FACT[q])).T@w
tot2 = tot - zc; print("minus pathC coincid", r(tot2-Dmc)); tot3 = tot2 - za; print("minus pathA coincid", r(tot3-Dmc))
