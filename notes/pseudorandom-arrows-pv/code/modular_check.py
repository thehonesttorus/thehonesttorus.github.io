# Modular characterization: u_pq (x) 1 is an eigenoperator of sigma_t = Ad(rho^{it}) iff nu_p = nu_q.
import numpy as np, itertools
N = 3
# correlated two-bit head, fair independent third bit:  nu(x1,x2) = [0.4, 0.1, 0.1, 0.4]
head = {(0,0): .4, (0,1): .1, (1,0): .1, (1,1): .4}
w = np.array([head[(x[0], x[1])] * 0.5 for x in itertools.product([0,1], repeat=N)])
rho = np.diag(w * 2**N)                       # density w.r.t. normalized trace
def unit(p, q, n):                            # u_pq (x) 1 on N bits
    P = np.zeros((2**n, 2**n)); P[int(''.join(map(str,p)),2), int(''.join(map(str,q)),2)] = 1
    return np.kron(P, np.eye(2**(N-n)))
t = 0.7
U = np.diag(np.exp(1j * t * np.log(np.diag(rho))))
for n in [1, 2]:
    for p in itertools.product([0,1], repeat=n):
        for q in itertools.product([0,1], repeat=n):
            if p >= q: continue
            u = unit(p, q, n); s = U @ u @ U.conj().T
            nz = np.abs(u) > 0
            ratios = s[nz] / u[nz]
            eig = np.allclose(ratios, ratios[0])
            # continuation laws
            def cont(p):
                idx = [i for i, x in enumerate(itertools.product([0,1], repeat=N)) if x[:n] == p]
                v = w[idx]; return v / v.sum()
            same = np.allclose(cont(p), cont(q))
            print(f"level {n}: arrow {p}->{q}: modular eigenoperator={eig}, continuation laws equal={same}"
                  + (f", eigenvalue^(1/it) = {np.exp(np.log(ratios[0]).imag / t):.3f} vs nu[p]/nu[q] = {w[[i for i,x in enumerate(itertools.product([0,1],repeat=N)) if x[:n]==p]].sum()/w[[i for i,x in enumerate(itertools.product([0,1],repeat=N)) if x[:n]==q]].sum():.3f}" if eig else ""))
