# Lumped Ornstein-Uhlenbeck kernel on the activation tiles of a random bias-free ReLU net:
#   K_s(t,t') = P(X_s in t' | X_0 in t),  X_s = e^{-s} X_0 + sqrt(1-e^{-2s}) Z.
# Theorem: every eigenvalue of K_s on mean-zero tile functions lies in [0, e^{-s}], for ANY partition.
import numpy as np
rng = np.random.default_rng(3)
def patterns(Ws, X):
    H = X; bits = []
    for W in Ws:
        Z = H @ W.T; bits.append(Z > 0); H = np.maximum(Z, 0)
    B = np.concatenate(bits, 1)
    return B
def tile_ids(B, keep=None):
    keys = np.packbits(B, axis=1)
    keys = keys.view(np.dtype((np.void, keys.shape[1]))).ravel()
    u, inv = np.unique(keys, return_inverse=True)
    return u, inv
for (d, n, L) in [(3, 6, 3), (4, 8, 3), (5, 10, 4)]:
    Ws = [rng.standard_normal((n, d)) * np.sqrt(2/d)] + [rng.standard_normal((n, n)) * np.sqrt(2/n) for _ in range(L-1)]
    N = 3_000_000
    for s in [0.1, 0.5, 1.0]:
        r = np.exp(-s)
        X0 = rng.standard_normal((N, d)); Xs = r*X0 + np.sqrt(1-r*r)*rng.standard_normal((N, d))
        B0, Bs = patterns(Ws, X0), patterns(Ws, Xs)
        u, inv = tile_ids(np.concatenate([B0, Bs]))
        a, b = inv[:N], inv[N:]
        m = np.bincount(np.concatenate([a, b]), minlength=len(u)) / (2*N)
        # merge rare tiles into one class (still a partition, so the theorem applies)
        big = m >= 2e-4
        lab = np.where(big, np.cumsum(big) - 1, big.sum())
        k = big.sum() + 1
        a, b = lab[a], lab[b]
        C = np.zeros((k, k)); np.add.at(C, (a, b), 1.0); C = (C + C.T) / (2*N)
        p = C.sum(1)
        M = C / np.sqrt(np.outer(p, p))
        ev = np.sort(np.linalg.eigvalsh(M))[::-1]
        print(f"d={d} n={n} L={L} s={s}: tiles={len(u)} (kept {k-1}+1), lambda_2={ev[1]:.4f} <= e^-s={r:.4f}, min eig={ev[-1]:.4f}")
