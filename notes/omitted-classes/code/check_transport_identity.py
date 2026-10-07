import numpy as np, itertools, math
src = open("yclasses.py").read()
pre_fn = src[src.index("def offd"):src.index("def truth")]
Wcol = np.load("../official/W_off9.npy").astype(np.float64); L, n, _ = Wcol.shape
ia = np.array([i for i in range(n) for j in range(n) if i != j]); ib = np.array([j for i in range(n) for j in range(n) if i != j])
exec(pre_fn)
# regenerate the samples of mcstats (N=400000, seed 11, one batch per half)
N = 400000; rng = np.random.default_rng(11); Wf = Wcol.astype(np.float32)
bs = []; done = 0
while done < N // 2:
    b = min(32768, N // 2 - done); bs.append(rng.standard_normal((n, b), dtype=np.float32)); done += b
Ys = [np.concatenate(bs, 1)]
F = np.load("mt_off9_h0.npz")
def k4tensor(X):
    X = X.astype(np.float64); X = X - X.mean(1, keepdims=True); m = X.shape[1]
    C = X @ X.T / m
    M4 = np.einsum("ap,bp,cp,dp->abcd", X, X, X, X, optimize=True) / m
    return M4 - np.einsum("ab,cd->abcd", C, C) - np.einsum("ac,bd->abcd", C, C) - np.einsum("ad,bc->abcd", C, C)
Y = Ys[0]; ys = []; zs = []
for l in range(L):
    Z = Wf[l] @ Y; zs.append(Z); Y = np.maximum(Z, 0); ys.append(Y)
for l in range(L - 1):
    W = Wcol[l + 1]; H = W * W
    Ky = k4tensor(ys[l]); Kz = k4tensor(zs[l + 1])
    # literal transport identity
    Kt = np.einsum("ia,jb,kc,ld,abcd->ijkl", W, W, W, W, Ky, optimize=True)
    # class masks of y tensor
    idx = np.array(list(itertools.product(range(n), repeat=4))).reshape(n, n, n, n, 4)
    nd = np.apply_along_axis(lambda v: len(set(v)), 4, idx)
    Krem = np.where(nd >= 3, Ky, 0.0)
    Rbf = np.einsum("ia,jb,kc,ld,abcd->ijkl", W, W, W, W, Krem, optimize=True)
    # slices from the npz vs brute force
    d_y = np.array([Ky[a, a, a, a] for a in range(n)]); K_y = np.array([[Ky[a, a, b, b] for b in range(n)] for a in range(n)]); np.fill_diagonal(K_y, 0)
    B_y = np.array([[Ky[a, a, a, b] for b in range(n)] for a in range(n)]); np.fill_diagonal(B_y, 0)
    e_npz = max(abs(F["k4_y"][l] - d_y).max(), abs(offd(F["K22_y"][l]) - K_y).max(), abs(offd(F["K31_y"][l]) - B_y).max())
    tz = (np.array([Kz[i, i, i, i] for i in range(n)]), np.array([Kz[i, i, j, j] for i, j in zip(ia, ib)]),
          offd(np.array([[Kz[i, i, i, j] for j in range(n)] for i in range(n)])))
    Tp = T_pair(W, H, d_y, K_y, B_y)
    R = tuple(a - b for a, b in zip(tz, Tp))
    Rb = (np.array([Rbf[i, i, i, i] for i in range(n)]), np.array([Rbf[i, i, j, j] for i, j in zip(ia, ib)]),
          offd(np.array([[Rbf[i, i, i, j] for j in range(n)] for i in range(n)])))
    print(l, "npz-vs-bf %.1e" % e_npz, "z=W#y %.1e" % (abs(Kt - Kz).max() / abs(Kz).max()),
          "R vs bf", ["%.1e" % (np.abs(a - b).max() / np.abs(b).max()) for a, b in zip(R, Rb)],
          "|R|/|t|", ["%.2f" % (np.linalg.norm(a) / np.linalg.norm(b)) for a, b in zip(Rb, tz)])
