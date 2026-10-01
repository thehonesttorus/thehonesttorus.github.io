"""Anatomy of the exact old slice (age > A) at width n: scale-mode share and the residual's modes."""
import sys, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "../../bench"))
import bench, costate
name, i, A = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
S = bench.load_set(name); W = bench.weights(S, i).astype(float)
rec = []
costate.predict(W, A=A, oldfilter=lambda X: X, record=rec, keep_old=True)
cos = lambda a, b: abs(a @ b) / np.linalg.norm(a) / np.linalg.norm(b)
for x in rec:
    if x.get("Sold") is None or x["layer"] % 2: continue
    So, K, v = x["Sold"], x["Kg"], x["vecs"]
    g = np.sum(K * So) / np.sum(K * K)
    R = So - g * K
    off = lambda X: X - np.diag(np.diag(X))
    sh = 1 - np.sum(R ** 2) / np.sum(So ** 2)
    sho = 1 - np.sum(off(R) ** 2) / np.sum(off(So) ** 2)
    Uu, s, Vt = np.linalg.svd(off(R)); e = s ** 2 / np.sum(s ** 2)
    al = " ".join(f"{k}:{cos(Uu[:, 0], w):.2f}/{cos(Vt[0], w):.2f}" for k, w in v.items())
    print(f"L{x['layer']:2d} gamma {g:.4f} | scale-mode share of |S_old|^2: all {sh:.2f} off {sho:.2f} | residual off top1 {e[0]:.2f} top4 {e[:4].sum():.2f} top16 {e[:16].sum():.2f} | res u1/v1: {al}", flush=True)
