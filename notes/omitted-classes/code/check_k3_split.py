# brute-force check of the kappa3 class split: R3 = transport of the all-distinct y class
import numpy as np, itertools
Wcol = np.load("../official/W_off9.npy").astype(np.float64); L, n, _ = Wcol.shape
N = 400000; rng = np.random.default_rng(11); Wf = Wcol.astype(np.float32)
bs = []; done = 0
while done < N // 2:
    b = min(32768, N // 2 - done); bs.append(rng.standard_normal((n, b), dtype=np.float32)); done += b
Y = np.concatenate(bs, 1); ys = []
for l in range(L):
    Z = Wf[l] @ Y; Y = np.maximum(Z, 0); ys.append(Y)
F = np.load("mt_off9_h0.npz")
for l in range(L - 1):
    W = Wcol[l + 1]; H = W * W; X = ys[l].astype(np.float64); X = X - X.mean(1, keepdims=True)
    K3 = np.einsum("ap,bp,cp->abc", X, X, X) / X.shape[1]
    nd = np.array([[[len({a, b, c}) for c in range(n)] for b in range(n)] for a in range(n)])
    Rt = np.einsum("ia,jb,kc,abc->ijk", W, W, W, np.where(nd == 3, K3, 0.0))
    R3bf = (np.array([Rt[i, i, i] for i in range(n)]), np.array([[Rt[i, i, c] if i != c else 0 for c in range(n)] for i in range(n)]))
    k3v = F["k3_y"][l].astype(np.float64); D = F["D21_y"][l].astype(np.float64).copy(); np.fill_diagonal(D, 0)
    G3 = D @ W.T
    d3 = (W * H) @ k3v + 3 * np.einsum("ia,ai->i", H, G3); d21 = (H * k3v[None, :]) @ W.T + H @ D @ W.T + 2 * (W * G3.T) @ W.T; np.fill_diagonal(d21, 0)
    t3 = (F["k3"][l + 1].astype(np.float64), F["D21"][l + 1].astype(np.float64).copy()); np.fill_diagonal(t3[1], 0)
    R3 = (t3[0] - d3, t3[1] - d21)
    print(l, ["%.1e" % (np.abs(a - b).max() / np.abs(b).max()) for a, b in zip(R3, R3bf)])
