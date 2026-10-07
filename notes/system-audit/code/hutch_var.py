# Stochastic diagonal estimation of the kappa4 pair-class residual operators, measured on the chain's stored states:
# for an operator R (n x n), one Rademacher probe Z = eps * (R eps) has E Z = diag R and Var Z_i = sum_{j != i} R_ij^2.
# Signal = the diagonal we want (the quenched part); noise = sqrt(mean_i Var Z_i); probes for unit SNR = (noise/signal)^2.
import sys, pickle, numpy as np
net = int(sys.argv[1]); srcs = [int(v) for v in sys.argv[2].split(",")]
Wc = np.load(f"../official/W_off{net}.npy").astype(np.float64); n = 1024
D = {d["layer"]: d for d in pickle.load(open(f"v29k4_off{net}.pkl", "rb"))["dumps"]}
for s in srcs:
    W = Wc[s + 1]; W2 = W * W
    K = D[s]["K22"].astype(np.float64).copy(); np.fill_diagonal(K, 0.0); K = 0.5 * (K + K.T)
    M = 3.0 * W2 @ K @ W2.T                                   # pair-class operator, diag = pair class
    m = W2.mean(1)                                            # row means: W o W = m 1^T + xi
    xi = W2 - m[:, None]
    Q = 3.0 * xi @ K @ xi.T                                   # deflated: the quenched operator, diag = quenched part
    lam, V = np.linalg.eigh(K); o = np.argsort(-np.abs(lam))[:4]
    K4 = (V[:, o] * lam[o]) @ V[:, o].T
    R4 = 3.0 * xi @ (K - K4) @ xi.T                           # remainder after the exact rank-4 part
    for tag, A in (("pair-class operator", M), ("deflated by row means", Q), ("remainder after rank 4", R4)):
        d = np.diag(A); off = A - np.diag(d)
        sig = np.sqrt(np.mean(d * d)); noise = np.sqrt(np.mean((off * off).sum(1)))
        print(f"source {s:2d} {tag:24s}: diag rms {sig:.3e}  probe noise rms {noise:.3e}  probes for SNR 1: {(noise / sig) ** 2:9.0f}"
              f"  (cost at 6 n^2 per probe: {6 * (noise / sig) ** 2 / 2048:.1f} units per layer)", flush=True)
