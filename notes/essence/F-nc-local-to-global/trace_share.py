"""Team F test T1: share of the exact old (2,1) slice (age > A) at width n that lies in the trace channel
{u x^T + y u^T} (u = scale vector), i.e. the defect sector of the fresh-weight commuting square.
Reuses costate's exact first-order co-state (all pairs) and its recorded old slices."""
import sys, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
CS = os.path.join(HERE, "../../fresh-slate/breakthrough/costate")
sys.path.insert(0, CS); sys.path.insert(0, os.path.join(CS, "../../bench"))
import bench, costate
name, i, A = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
S = bench.load_set(name); W = bench.weights(S, i).astype(float)
rec = []
costate.predict(W, A=A, oldfilter=lambda X: X, record=rec, keep_old=True)
off = lambda X: X - np.diag(np.diag(X))
rng = np.random.default_rng(0)
def chan(X, u):
    u = u / np.linalg.norm(u); a = u @ X; b = X @ u
    return np.outer(u, a) + np.outer(b, u) - (u @ X @ u) * np.outer(u, u)
for x in rec:
    if x.get("Sold") is None: continue
    So = off(x["Sold"]); v = x["vecs"]; n = So.shape[0]; E = np.sum(So ** 2)
    K = off(x["Kg"]); g = np.sum(K * So) / np.sum(K * K); scale = 1 - np.sum((So - g * K) ** 2) / E
    out = [f"L{x['layer']:2d} scale-mode {scale:.3f}"]
    for k in ("s", "s2", "one", "m", "Phi"):
        P = chan(So, v[k]); out.append(f"trace[{k}] {np.sum(P**2)/E:.3f}")
    # two-vector channel: u in span{s2, m}
    Q, _ = np.linalg.qr(np.stack([v["s2"], v["m"]], 1))
    P2 = Q @ (Q.T @ So) + (So @ Q) @ Q.T - Q @ (Q.T @ So @ Q) @ Q.T
    out.append(f"trace[s2,m] {np.sum(P2**2)/E:.3f}"); out.append(f"E_old {E:.4e} E_trace {np.sum(P2**2):.4e} E_bulk {E-np.sum(P2**2):.4e} E_young {np.sum(off(x['Syoung'])**2) if x.get('Syoung') is not None else float('nan'):.4e}")
    R2 = So - P2; Kr = K - (Q @ (Q.T @ K) + (K @ Q) @ Q.T - Q @ (Q.T @ K @ Q) @ Q.T)
    gu = np.sum(Kr * R2) / np.sum(Kr * Kr); out.append(f"union[scale+trace] {1 - np.sum((R2 - gu * Kr) ** 2) / E:.3f}")
    r = rng.standard_normal(n); out.append(f"control[rand] {np.sum(chan(So, r)**2)/E:.4f}")
    print(" | ".join(out), flush=True)
