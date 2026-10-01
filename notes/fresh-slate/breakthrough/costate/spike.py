"""Anatomy of the old (age > A) slice at each layer: top singular vectors vs reference-chain vectors."""
import sys, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "../../bench"))
sys.path.insert(0, os.path.join(HERE, "../../designs/heisenberg"))
import bench, costate
from hd import Gauss
name, i, A = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 1
S = bench.load_set(name); W = bench.weights(S, i).astype(float)
recA, recF = [], []
costate.predict(W, A=A, record=recA); costate.predict(W, record=recF)
# reference-chain vectors per layer (closure chain is close enough for directions)
L, n, _ = W.shape
m = np.zeros(n); C = W[0].T @ W[0]; vecs = []
for l in range(L):
    G = Gauss(m, C, K=2)
    vecs.append(dict(one=np.ones(n), m=m.copy(), Ea=G.Ea.copy(), Phi=G.Phi.copy(), s=G.s.copy(), p0=G.Ed[0].copy()))
    if l + 1 < L:
        Wn = W[l + 1]; m = G.Ea @ Wn; C = Wn.T @ G.cov_a() @ Wn
cos = lambda a, b: abs(a @ b) / np.linalg.norm(a) / np.linalg.norm(b)
for l in range(2, L, 2):
    So = recF[l]["S"] - recA[l]["S"]
    off = So - np.diag(np.diag(So))
    Uu, s, Vt = np.linalg.svd(off)
    e = s ** 2 / np.sum(s ** 2)
    v = vecs[l]
    best = {k: (cos(Uu[:, 0], x), cos(Vt[0], x)) for k, x in v.items() if np.linalg.norm(x) > 0}
    bs = " ".join(f"{k}:{a:.2f}/{b:.2f}" for k, (a, b) in best.items())
    print(f"layer {l:2d} |off|/|diag| {np.linalg.norm(off)/np.linalg.norm(np.diag(So)):.1f}  energy top1 {e[0]:.2f} top8 {e[:8].sum():.2f}  cos(u1,v1 vs x): {bs}")
