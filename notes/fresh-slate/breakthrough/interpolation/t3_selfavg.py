"""T3: how much of the closure error e_j (final layer, n=1024) is a smooth function of per-neuron quenched features
(alpha_j = m_j/s_j, s_j)?  = the self-averaging (chaos-0, disorder-conditional) part.  Leave-one-MLP-out fit.
Also: closure-predicted Var(|a-mu|^2)/(E|a-mu|^2)^2 = 2||C||_F^2/tr(C)^2 per layer (compare t1's measured Var tau)."""
import numpy as np, sys
from common import bench, phi, Phi, chi_mean_ratio
from t2_closures import hk

def closure_feats(W, K=2):
    W = W.astype(np.float64); L, n, _ = W.shape
    m = np.zeros(n); S = W[0].T @ W[0]; vt = []
    for l in range(L):
        if l > 0: m = mu @ W[l]; S = W[l].T @ C @ W[l]
        v = np.diag(S); s = np.sqrt(v); a = m / s; P = Phi(a); p = phi(a)
        mu = m * P + s * p; sec = (m * m + v) * P + m * s * p
        R = S / np.outer(s, s); H = hk(a, K)
        C = np.outer(s * H[0], s * H[0]) * R + 0.5 * np.outer(s * H[1], s * H[1]) * R * R
        np.fill_diagonal(C, sec - mu * mu)
        vt.append(2 * (C ** 2).sum() / np.trace(C) ** 2)
    return mu * chi_mean_ratio(n), a, s, np.array(vt)

S = bench.load_set("w1024_d16"); n = S["width"]
X, Y, G = [], [], []
for i in range(len(S["seeds"])):
    W = bench.weights(S, i); mu, a, s, vt = closure_feats(W)
    if i == 0: print("closure 2||C||^2/tr^2 in units 2/n:", np.round(vt * n / 2, 2))
    e = S["means"][i][-1] - mu
    X.append(np.stack([a, s], 1)); Y.append(e); G.append(np.full(n, i))
X = np.concatenate(X); Y = np.concatenate(Y); G = np.concatenate(G)
def basis(X, deg):
    a, s = X[:, 0], X[:, 1]
    cols = [s * a ** k for k in range(deg + 1)] + [s * phi(a) * a ** k for k in range(deg + 1)] + [a ** k for k in range(deg + 1)]
    return np.stack(cols, 1)
for deg in (0, 1, 2, 4):
    res = np.empty_like(Y)
    for g in np.unique(G):
        tr = G != g; B = basis(X, deg)
        c, *_ = np.linalg.lstsq(B[tr], Y[tr], rcond=None)
        res[~tr] = Y[~tr] - B[~tr] @ c
    print(f"deg {deg}: MSE before {np.mean(Y**2):.3e}  after (LOO-MLP) {np.mean(res**2):.3e}  "
          + " ".join(f"{np.mean(res[G==g]**2):.2e}" for g in np.unique(G)))
