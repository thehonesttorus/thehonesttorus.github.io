"""Team F test T1/T2 at the readout: inject only parts of the exact old (2,1) slice (age > A) into costate's
first-order co-state at n = 1024 and score raw final MSE. Filters: exact diagonal plus
  cap   : off-diagonal projected on the cap channel trace[s2,m] united with the scale mode;
  capR k: cap + rank-k truncation of the remaining bulk;
  rankR k (for comparison): plain rank-k of the whole off-diagonal (costate's filter).
Paired against 'all' (identity) and 'none'."""
import sys, os, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
CS = os.path.join(HERE, "../../fresh-slate/breakthrough/costate")
sys.path.insert(0, CS); sys.path.insert(0, os.path.join(CS, "../../bench"))
import bench, costate
name, i, A = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]); which = sys.argv[4].split(",")
S = bench.load_set(name); W = bench.weights(S, i).astype(float); truth = S["means"][i]; noise = S["noise"][i]
offd = lambda X: X - np.diag(np.diag(X))
def chanproj(X, Q): return Q @ (Q.T @ X) + (X @ Q) @ Q.T - Q @ (Q.T @ X @ Q) @ Q.T
def make(kind, k, rec):
    def f(So):
        x = rec[-1]; v = x["vecs"]; O = offd(So)
        if kind == "none": return np.diag(np.diag(So))
        if kind == "all": return So
        if kind == "rankR":
            U, s, Vt = np.linalg.svd(O); return np.diag(np.diag(So)) + offd((U[:, :k] * s[:k]) @ Vt[:k])
        Q, _ = np.linalg.qr(np.stack([v["s2"], v["m"]], 1))
        P = chanproj(O, Q); B = O - P; K = offd(x["Kg"]); K = K - chanproj(K, Q)
        g = np.sum(K * B) / np.sum(K * K); cap = P + g * K; B = B - g * K
        out = np.diag(np.diag(So)) + cap
        if kind == "capR" and k > 0:
            U, s, Vt = np.linalg.svd(B); out = out + offd((U[:, :k] * s[:k]) @ Vt[:k])
        return out
    return f
for w in which:
    kind, k = (w.split(":") + ["0"])[:2]; k = int(k); rec = []
    pred = costate.predict(W, A=A, oldfilter=make(kind, k, rec), record=rec, keep_old=True)
    raw = float(((pred[-1] - truth[-1]) ** 2).mean() - noise)
    print(json.dumps(dict(set=name, mlp=i, A=A, filter=w, raw=raw)), flush=True)
