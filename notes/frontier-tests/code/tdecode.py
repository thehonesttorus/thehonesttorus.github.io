# Transported decode: accumulate candidate local covariance defects M_k(l') through the actual layers
# (linear response: offdiag F1F1^T o (W A W^T), diag dVar/dvar * diag(W A W^T)) and regress the one-step residual
# r_l on sens_i * w_i^T A_k(l-1) w_i, pooled over layers 5..15. Fit on one net, evaluate on another.
import pickle, sys, numpy as np
from scipy.stats import norm
def load(net):
    D = {d["layer"]: d for d in pickle.load(open(f"dump_oracleall_off{net}.pkl", "rb"))}
    pred = np.load(f"pred_oracleall_off{net}.npy"); mt = np.load(f"../official/truth_off{net}.npz")["m"].astype(np.float64)
    Wc = np.load(f"../official/W_off{net}.npy").astype(np.float32)
    return D, pred, mt, Wc
def post(d):
    s = np.sqrt(d["var"]); a = d["mu"] / s; ph, Ph = norm.pdf(a), norm.cdf(a)
    m = s * (a * Ph + ph); v = s**2 * ((a**2 + 1) * Ph + a * ph) - m**2
    # dVar(relu)/dvar at fixed mean: dE[relu^2]/dvar = Phi, dm/dvar = phi/(2s)
    dv = Ph - 2 * m * ph / (2 * s)
    return s, a, ph, Ph, m, v, dv
def offd(M): M = M.copy(); np.fill_diagonal(M, 0.0); return M
NAMES = ["lin", "mehl2", "g4rank1", "k4iiij", "d21", "diagv", "mm"]
def local_defects(d):
    s, a, ph, Ph, m, v, dv = post(d)
    C = offd(d["C_pre"]); F1, F2, F3 = Ph, ph / s, -a * ph / s**2
    out = {"lin": np.outer(F1, F1) * C, "mehl2": 0.5 * np.outer(F2, F2) * C**2,
           "g4rank1": offd(np.outer(F2 * s**2, F2 * s**2)), "k4iiij": (np.outer(F1, F3 * s**2) + np.outer(F3 * s**2, F1)) * C,
           "d21": offd(np.outer(F1, F2) * d["D21"].T + np.outer(F2, F1) * d["D21"]) if d["D21"] is not None else np.zeros_like(C),
           "diagv": np.diag(v), "mm": np.outer(m, m)}
    return {k: v_.astype(np.float32) for k, v_ in out.items()}
def features(net, l0=5):
    D, pred, mt, Wc = load(net)
    acc = None; X, Y = [], []
    layers = sorted(D)
    for l in layers:
        if l == 15 or D[l]["C_pre"] is None:
            pass
        # features for layer l use acc at l-1 (post)
        if acc is not None and l >= l0:
            W = Wc[l]; c = D[l]
            sl = np.sqrt(c["var"]); al = c["mu"] / sl; sens = norm.pdf(al) / (2 * sl)
            F = np.column_stack([sens * np.einsum("ij,ij->i", W @ acc[k], W) for k in NAMES])
            X.append(F); Y.append(mt[l] - pred[l])
        if D[l]["C_pre"] is None: break
        # transport acc from l-1 post to l post, then add local defects at l
        loc = local_defects(D[l])
        if acc is None:
            acc = loc
        else:
            W = Wc[l]; s, a, ph, Ph, m, v, dv = post(D[l])
            new = {}
            for k in NAMES:
                Q = W @ acc[k] @ W.T
                dq = np.diag(Q).copy()
                T = np.outer(Ph, Ph).astype(np.float32) * Q
                np.fill_diagonal(T, dv * dq)
                new[k] = T + loc[k]
            acc = new
    return np.vstack(X), np.concatenate(Y), len(X)
Xa, Ya, La = features(int(sys.argv[1])); Xb, Yb, Lb = features(int(sys.argv[2]))
def fit(X, Y):
    X1 = np.column_stack([np.ones(len(Y)), X]); b, *_ = np.linalg.lstsq(X1, Y, rcond=None); return b
def r2(b, X, Y):
    X1 = np.column_stack([np.ones(len(Y)), X]); return 1 - np.mean((Y - X1 @ b)**2) / np.mean(Y**2)
b = fit(Xa, Ya)
print("pooled layers 5..15; R2 = 1 - MSE(resid)/MSE(r)")
print(f"joint: in-sample (net {sys.argv[1]}) {r2(b, Xa, Ya):.3f}   out-of-sample (net {sys.argv[2]}) {r2(b, Xb, Yb):.3f}")
for j, k in enumerate(NAMES):
    bj = fit(Xa[:, [j]], Ya); print(f"  {k:8s} in {r2(bj, Xa[:, [j]], Ya):.3f}  out {r2(bj, Xb[:, [j]], Yb):.3f}  coef {bj[1]:+.3e}")
