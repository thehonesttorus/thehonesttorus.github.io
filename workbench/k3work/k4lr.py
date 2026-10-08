# How much of the quenched kappa4 pair-class term 3 diag(W2 K22 W2^T) - (mean field) does a rank-k approximation of
# the post-activation (2,2) slice K22 carry?  Exact eigenpairs and a randomized range finder (the chain's recipe:
# Gaussian sketch, one power iteration), on the production chain's dumped K22 at the given layers.
import sys, pickle, numpy as np
net = int(sys.argv[1]); srcs = [int(v) for v in sys.argv[2].split(",")]
Wc = np.load(f"../official/W_off{net}.npy").astype(np.float64); n = 1024
D = {d["layer"]: d for d in pickle.load(open(f"v29k4_off{net}.pkl", "rb"))["dumps"]}
cA, cI = 6.0 / (n + 4.0), -3.0 / ((n + 2.0) * (n + 4.0))
def mf(W2, K4, K):
    Kz = K - np.diag(np.diag(K))
    return 2.0 * W2 @ (cA * (K4 + Kz.sum(1)) + cI * (K4.sum() + Kz.sum()))
def exact_lr(W2, W4, K4, lam, V):
    Kd = (V * V) @ lam                                    # diagonal of the low-rank K, removed (pair class is a != b)
    return W4 @ K4 + 3.0 * ((W2 @ V) ** 2) @ lam - 3.0 * W4 @ Kd
rng = np.random.default_rng(0)
for s in srcs:
    W = Wc[s + 1]; W2 = W * W; W4 = W2 * W2
    K = D[s]["K22"].astype(np.float64).copy(); np.fill_diagonal(K, 0.0); K = 0.5 * (K + K.T); K4 = D[s]["K4v"].astype(np.float64)
    ex = W4 @ K4 + 3.0 * np.einsum("ij,ij->i", W2 @ K, W2); q = ex - mf(W2, K4, K)
    lam, V = np.linalg.eigh(K); o = np.argsort(-np.abs(lam)); lam, V = lam[o], V[:, o]
    en = np.cumsum(lam ** 2) / np.sum(lam ** 2)
    print(f"source {s}: |K22|_F {np.linalg.norm(K):.4f}; quenched rms {np.sqrt(np.mean(q*q)):.2e}; spectral energy in top k: " +
          " ".join(f"{k}:{en[k-1]:.3f}" for k in (1, 2, 4, 8, 16, 32, 64, 128, 256)), flush=True)
    for k in (1, 2, 4, 8, 16, 32, 64, 128):
        lk, Vk = lam[:k], V[:, :k]
        Kk = (Vk * lk) @ Vk.T
        qk = exact_lr(W2, W4, K4, lk, Vk) - mf(W2, K4, Kk)
        # randomized range finder of rank k (+8 oversampling), one power iteration
        Om = rng.standard_normal((n, k + 8)); Y = K @ Om; Qr, _ = np.linalg.qr(Y); Y = K @ (K @ Qr); Qr, _ = np.linalg.qr(Y)
        Bm = Qr.T @ K @ Qr; lb, Ub = np.linalg.eigh(Bm); ob = np.argsort(-np.abs(lb))[:k]; lb, Vb = lb[ob], Qr @ Ub[:, ob]
        qr_ = exact_lr(W2, W4, K4, lb, Vb) - mf(W2, K4, (Vb * lb) @ Vb.T)
        f = lambda a: (np.corrcoef(a, q)[0, 1], np.sqrt(np.mean((a - q) ** 2)) / np.sqrt(np.mean(q * q)))
        c1, e1 = f(qk); c2, e2 = f(qr_)
        print(f"   k {k:3d}: exact eigenpairs corr {c1:.3f} rel err {e1:.3f} | randomized corr {c2:.3f} rel err {e2:.3f} "
              f"| cost ~{3 * (k + 8) / 1024 * 2:.2f} units", flush=True)
