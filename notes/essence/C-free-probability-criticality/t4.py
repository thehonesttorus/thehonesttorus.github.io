"""T4: is the sub-free PR of old propagators the rank-one Phi-W (Perron) coupling?
Propagator U_{0->k} = W_1 diag(g_1) W_2 ... with g = mean gates (Gaussian-closure chain), the same gates randomly
permuted per layer (kills the coupling to W; free/traffic-independent diagonal), and pathwise 0/1 gates of one input.
python t4.py <set> <mlp>"""
import sys, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "../../fresh-slate/bench")); sys.path.insert(0, os.path.join(HERE, "../../fresh-slate/designs/heisenberg"))
import bench
from hd import Gauss
name, i = sys.argv[1], int(sys.argv[2])
S = bench.load_set(name); Ws = bench.weights(S, i).astype(np.float64); L, n, _ = Ws.shape
m = np.zeros(n); C = Ws[0].T @ Ws[0]; Phis = []
for l in range(L - 1):
    G = Gauss(m, C, K=8); Phis.append(G.Phi.copy())
    m = G.Ea @ Ws[l + 1]; C = Ws[l + 1].T @ G.cov_a() @ Ws[l + 1]
rng = np.random.default_rng(0)
x = rng.standard_normal(n); z = x @ Ws[0]; paths = []
for l in range(L - 1):
    paths.append((z > 0).astype(float)); z = np.maximum(z, 0) @ Ws[l + 1]
def pr(U):
    G = U.T @ U; t = np.trace(G); return n * t * t / np.sum(G * G) / n, t / n
s0 = 1
for lab, gates in (("mean gates", Phis), ("permuted mean gates", [rng.permutation(p) for p in Phis]), ("pathwise gates (1 input)", paths)):
    U = Ws[s0 + 0].copy(); out = []
    for k in range(s0, L - 1):
        a = k - s0
        p, t = pr(U)
        r = np.mean(gates[k] ** 4) / np.mean(gates[k] ** 2) ** 2
        out.append(f"{a}:{p * 2 * (a + 1):.3f}")
        U = (U * gates[k][None, :]) @ Ws[k + 1]
    print(f"{lab:26s} PR*2(age+1)/n by age: " + " ".join(out), flush=True)
print("r(Phi^2) by layer:", " ".join(f"{np.mean(p**4)/np.mean(p**2)**2:.2f}" for p in Phis))
