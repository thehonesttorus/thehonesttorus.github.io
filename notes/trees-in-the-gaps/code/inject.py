# Fresh gain injection of one layer from an exactly Gaussian source: MC versus the O(n^2) formula.
import numpy as np
from closure import closure
from residue import closure_with_residue
n = 256; Ws = list(np.load("W_n256_L16_s0.npy")); oc = closure(Ws, keep=True)
comp = {t: np.array([x["gamma"] for x in closure_with_residue(Ws, terms=(t,))]) for t in ["k4","F0","k3mu","Kpmu"]}
rng = np.random.default_rng(2)
def gamma_of(Z):
    mu = Z.mean(0); M2 = Z.T @ Z/len(Z); Ez2 = (Z*Z).mean(0); E22 = (Z*Z).T @ (Z*Z)/len(Z)
    k = E22 - np.outer(Ez2, Ez2) - 2*M2**2 + 2*np.outer(mu**2, mu**2)
    off = ~np.eye(Z.shape[1], dtype=bool); return np.mean((k/np.outer(Ez2, Ez2))[off])
for l in [2, 6, 10, 14]:
    d = oc[l]; Sig = d["R"]*np.outer(d["sig"], d["sig"]); Lc = np.linalg.cholesky(Sig + 1e-10*np.eye(n))
    gs = []
    for rep in range(4):
        y = d["mu"] + rng.standard_normal((100000, n)) @ Lc.T
        z = np.maximum(y, 0) @ Ws[l+1].T; gs.append(gamma_of(z))
    f = {t: comp[t][l+1]-comp[t][l] for t in comp}
    print(f"layer {l+1}->{l+2}: MC fresh gamma {np.mean(gs):+.5f} +- {np.std(gs)/2:.5f} | formula {sum(f.values()):+.5f} = " + " ".join(f"{t}:{f[t]:+.5f}" for t in f))
