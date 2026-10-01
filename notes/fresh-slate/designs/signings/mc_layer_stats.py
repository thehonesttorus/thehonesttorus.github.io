"""per-layer MC statistics of pre-activations (mean, var, k3, k4, covariance, (2,1) slice) for one baked MLP."""
import sys, numpy as np
w, seed, N = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3])
W = np.load(f"/root/sg/truth/W_w{w}_s{seed}.npy").astype(np.float64); L, n, _ = W.shape
d = np.load(f"/root/sg/truth/w{w}_s{seed}.npz"); amean = d["S"] / d["n"]
zmean = np.stack([np.zeros(n)] + [amean[l - 1] @ W[l] for l in range(1, L)])
rng = np.random.default_rng(77); S = np.zeros((L, 4, n)); C = np.zeros((L, n, n)); C21 = np.zeros((L, n, n)); done = 0
while done < N:
    x = rng.standard_normal((1 << 15, n)); a = x
    for l in range(L):
        z = a @ W[l]; a = np.maximum(z, 0); zc = z - zmean[l]
        z2 = zc * zc
        S[l, 0] += zc.sum(0); S[l, 1] += z2.sum(0); S[l, 2] += (z2 * zc).sum(0); S[l, 3] += (z2 * z2).sum(0)
        C[l] += zc.T @ zc; C21[l] += z2.T @ zc
    done += 1 << 15
m = S / done; mu = m[:, 0]
var = m[:, 1] - mu ** 2
k3 = m[:, 2] - 3 * mu * m[:, 1] + 2 * mu ** 3
k4 = m[:, 3] - 4 * mu * m[:, 2] + 6 * mu ** 2 * m[:, 1] - 3 * mu ** 4 - 3 * var ** 2
E2 = C / done; Cz = E2 - mu[:, :, None] * mu[:, None, :]
E21 = C21 / done
D = E21 - m[:, 1][:, :, None] * mu[:, None, :] - 2 * E2 * mu[:, None, :] + 2 * (mu ** 2)[:, :, None] * mu[:, None, :]
np.savez(f"/root/sg/layerstats_w{w}_s{seed}.npz", mz=zmean + mu, var=var, k3=k3, k4=k4, Cz=Cz, D=D, N=done)
print("done", done)
