"""T15: MC (bake-centred) of C_l = E[a~ a~^T], Q_l = E[|a~|^2 a~ a~^T], E|a~|^2, E|a~|^4 at every layer (post-activation),
-> M_l = Q_l - tr(C_l) C_l - 2 C_l^2 (kappa4 matrix trace), X_l = Var|a~|^2 - 2||C_l||^2.  Saves float32 to scratch."""
import sys, numpy as np, time, os
from common import bench
name, i, N = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]); chunk = 4096
out = sys.argv[4]
S = bench.load_set(name); W = bench.weights(S, i); T = S["means"][i].astype(np.float32); L, n, _ = W.shape
C = np.zeros((L, n, n)); Q = np.zeros((L, n, n)); r2 = np.zeros(L); r4 = np.zeros(L)
rng = np.random.default_rng(31); done = 0; t0 = time.time()
while done < N:
    h = rng.standard_normal((chunk, n)).astype(np.float32)
    for l in range(L):
        h = np.maximum(h @ W[l], 0.0); u = h - T[l]
        q = (u * u).sum(1)
        C[l] += u.T @ u; Q[l] += (u * q[:, None]).T @ u
        r2[l] += q.sum(dtype=np.float64); r4[l] += (q.astype(np.float64) ** 2).sum()
    done += chunk
    if done % (chunk * 16) == 0: print(done, f"{time.time()-t0:.0f}s", flush=True)
C /= done; Q /= done; r2 /= done; r4 /= done
M = np.stack([Q[l] - np.trace(C[l]) * C[l] - 2 * C[l] @ C[l] for l in range(L)])
X = np.array([r4[l] - r2[l] ** 2 - 2 * (C[l] ** 2).sum() for l in range(L)])
np.savez(out, C=C.astype(np.float32), M=M.astype(np.float32), X=X, N=done)
print("X:", X, flush=True)
