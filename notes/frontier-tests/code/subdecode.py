# Subspace decode: if the state error at l-1 is concentrated in a k-dim subspace U, then
# r_l(i) ~ sens_i * g_i^T Delta g_i (covariance channel) + skew_i * <g_i^{x3}, T> (kappa3 channel), g_i = U^T w_i.
# Regress r_l on the k(k+1)/2 quadratic (and cubic) features; compare with chance (p/n) and random subspaces.
import pickle, sys, itertools, numpy as np
from scipy.stats import norm
net = int(sys.argv[1]) if len(sys.argv) > 1 else 0
D = {d["layer"]: d for d in pickle.load(open(f"dump_oracleall_off{net}.pkl", "rb"))}
pred = np.load(f"pred_oracleall_off{net}.npy"); mt = np.load(f"../official/truth_off{net}.npz")["m"].astype(np.float64)
Wcol = np.load(f"../official/W_off{net}.npy").astype(np.float64)
rng = np.random.default_rng(1)
def r2(X, r):
    X = np.column_stack([np.ones(len(r)), X]); beta, *_ = np.linalg.lstsq(X, r, rcond=None)
    return 1 - np.var(r - X @ beta) / np.var(r)
def feats(G, w, deg):
    k = G.shape[1]; cols = []
    for a, b in itertools.combinations_with_replacement(range(k), 2):
        cols.append(w[0] * G[:, a] * G[:, b])
    if deg == 3:
        for a, b, c in itertools.combinations_with_replacement(range(k), 3):
            cols.append(w[1] * G[:, a] * G[:, b] * G[:, c])
    return np.column_stack(cols)
print("layer | k=4 quad: top / rand | k=8 quad: top / rand | k=6 quad+cub: top / rand   (chance = p/1024)")
for l in sorted(D):
    if l - 1 not in D or D[l - 1]["C_pre"] is None: continue
    p, c = D[l - 1], D[l]
    s = np.sqrt(p["var"]); a = p["mu"] / s; ph, Ph = norm.pdf(a), norm.cdf(a)
    m = s * (a * Ph + ph); vpost = s**2 * ((a**2 + 1) * Ph + a * ph) - m**2
    Coff = p["C_pre"].copy(); np.fill_diagonal(Coff, 0)
    Cpost = np.outer(Ph, Ph) * Coff + 0.5 * np.outer(ph / s, ph / s) * Coff**2 + np.diag(vpost)
    ev, U = np.linalg.eigh(Cpost); U = U[:, ::-1]
    W = Wcol[l]; r = mt[l] - pred[l]
    sl = np.sqrt(c["var"]); al = c["mu"] / sl
    w = (norm.pdf(al) / (2 * sl), -al * norm.pdf(al) / sl**2 / 6)
    out = []
    for k, deg in [(4, 2), (8, 2), (6, 3)]:
        Gt = W @ U[:, :k]
        Q, _ = np.linalg.qr(rng.standard_normal((1024, k))); Gr = W @ Q
        Xt, Xr = feats(Gt, w, deg), feats(Gr, w, deg)
        out.append(f"{r2(Xt, r):.3f} / {r2(Xr, r):.3f} (ch {Xt.shape[1]/1024:.3f})")
    print(f"{l:4d}  | " + " | ".join(out) + f"   top-8 eig share {ev[::-1][:8].sum()/ev.sum():.2f}")
