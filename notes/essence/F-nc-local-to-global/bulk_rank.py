"""Team F test for T2: is the through-string bulk of the old (2,1) slice (old slice minus the cap channel
trace[s2,m] and the scale mode) compressible? Top-k singular-value energy share of the bulk, against an
i.i.d. Gaussian n x n control and against the full old slice."""
import sys, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
CS = os.path.join(HERE, "../../fresh-slate/breakthrough/costate")
sys.path.insert(0, CS); sys.path.insert(0, os.path.join(CS, "../../bench"))
import bench, costate
name, i, A = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
S = bench.load_set(name); W = bench.weights(S, i).astype(float)
rec = []; costate.predict(W, A=A, oldfilter=lambda X: X, record=rec, keep_old=True)
off = lambda X: X - np.diag(np.diag(X))
ks = (1, 16, 64, 256, 512)
def share(M):
    s = np.linalg.svd(M, compute_uv=False) ** 2; c = np.cumsum(s) / s.sum(); return [c[k - 1] for k in ks]
n = W.shape[1]; ctrl = share(np.random.default_rng(0).standard_normal((n, n)))
print("k:", ks, " gaussian control:", " ".join(f"{x:.3f}" for x in ctrl))
for x in rec:
    if x.get("Sold") is None or x["layer"] not in (6, 10, 15): continue
    So = off(x["Sold"]); v = x["vecs"]
    Q, _ = np.linalg.qr(np.stack([v["s2"], v["m"]], 1))
    B = So - (Q @ (Q.T @ So) + (So @ Q) @ Q.T - Q @ (Q.T @ So @ Q) @ Q.T)
    K = off(x["Kg"]); K = K - (Q @ (Q.T @ K) + (K @ Q) @ Q.T - Q @ (Q.T @ K @ Q) @ Q.T)
    B = B - np.sum(K * B) / np.sum(K * K) * K
    sym = np.sum(((B + B.T) / 2) ** 2) / np.sum(B ** 2)
    print(f"L{x['layer']:2d} old:  " + " ".join(f"{x:.3f}" for x in share(So)) + f" | bulk ({np.sum(B**2)/np.sum(So**2):.3f} of old): " + " ".join(f"{x:.3f}" for x in share(B)) + f" | bulk sym share {sym:.3f}", flush=True)
