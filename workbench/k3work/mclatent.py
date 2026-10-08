# One Monte Carlo chunk for the latent-geometry test of the fourth-cumulant slices.
#   python mclatent.py NET NSAMP CHUNK K LAYERS OUTFILE      (LAYERS comma-separated, e.g. 3,6,9,11,13,15)
# At each listed layer l, with mu_l and the top-K eigenvectors U of the true covariance C_l (mc2_off{NET}_full.npz),
# every sample's pre-activation fluctuation d = z_l - mu_l is split into its collective part c = U U^T d and the
# residual r = d - c. For each of d, c, r it accumulates the raw-moment sums mcstats.py uses for the diagonal
# cumulants and the (1,1), (2,1), (3,1) slices: s1..s4 (per unit) and m11 = sum v v^T, m21 = sum v^2 v^T,
# m31 = sum v^3 v^T. The (3,1) slice of c is the latent model's prediction T4[U_i, U_i, U_i, U_j] with the true latent
# law (t = U^T d), so comparing the three variables splits each true slice into collective, residual and cross parts.
# Samples come from default_rng([4711, NET, CHUNK]); chunks are independent.
import sys, time, numpy as np
net, N, chunk, K = int(sys.argv[1]), int(float(sys.argv[2])), int(sys.argv[3]), int(sys.argv[4])
LAY = [int(x) for x in sys.argv[5].split(",")]; outf = sys.argv[6]
Wcol = np.load(f"../official/W_off{net}.npy").astype(np.float32)          # z_l = Wcol[l] @ y_(l-1)
L, n, _ = Wcol.shape
T = np.load(f"mc2_off{net}_full.npz")
MU = {l: np.asarray(T["mu"][l], np.float32) for l in LAY}
UK = {}
for l in LAY:
    ev, U = np.linalg.eigh(np.asarray(T["cov"][l], np.float64)); UK[l] = np.ascontiguousarray(U[:, ::-1][:, :K], np.float32)
VARS = ("d", "c", "r")
acc = {(l, v): dict(s1=np.zeros(n), s2=np.zeros(n), s3=np.zeros(n), s4=np.zeros(n),
                    m11=np.zeros((n, n)), m21=np.zeros((n, n)), m31=np.zeros((n, n))) for l in LAY for v in VARS}
rng = np.random.default_rng([4711, net, chunk])
B = 16384; done = 0; t0 = time.time()
while done < N:
    b = min(B, N - done)
    Y = rng.standard_normal((n, b), dtype=np.float32)
    for l in range(max(LAY) + 1):
        Z = Wcol[l] @ Y
        if l in MU:
            D = Z - MU[l][:, None]
            Cc = UK[l] @ (UK[l].T @ D); R = D - Cc
            for v, X in (("d", D), ("c", Cc), ("r", R)):
                a = acc[(l, v)]; X2 = X * X
                a["s1"] += X.sum(1, dtype=np.float64); a["s2"] += X2.sum(1, dtype=np.float64)
                a["s3"] += (X2 * X).sum(1, dtype=np.float64); a["s4"] += (X2 * X2).sum(1, dtype=np.float64)
                a["m11"] += X @ X.T; a["m21"] += X2 @ X.T; a["m31"] += (X2 * X) @ X.T
        Y = np.maximum(Z, 0.0, out=Z)
    done += b
    if (done // B) % 8 == 0:
        print(f"{done} samples, {time.time() - t0:.0f}s", flush=True)
np.savez(outf, c=done, layers=np.array(LAY), K=K,
         **{f"{k}_{v}_{l}": acc[(l, v)][k].astype(np.float64) for l in LAY for v in VARS for k in acc[(l, v)]})
print(f"done {done} samples in {time.time() - t0:.0f}s", flush=True)
