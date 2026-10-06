# The learner-independent bound on any address-based control variate: for an address T = P'x of dimension k, the
# variance a control variate c(T) can remove is at most Var(E[F|T]).  Measured on official net 0 (final-layer output
# F in R^n) by binning (k = 1, 2, 3; exact up to bin width) and by nearest-neighbour pairs (k = 8, 32; a lower bound
# on the conditional variance = upper bound on the explained fraction, which is the direction that favours the proposal).
# Addresses: top-k right singular vectors of W1 (the first layer's amplified directions) and random directions.
import numpy as np, sys, time
net = int(sys.argv[1]) if len(sys.argv) > 1 else 0; T = int(float(sys.argv[2])) if len(sys.argv) > 2 else 1_000_000
Wc = np.load(f"../official/W_off{net}.npy"); n, L = Wc.shape[1], Wc.shape[0]
Wf = [np.ascontiguousarray(W.T).astype(np.float32) for W in Wc]
_, _, Vt = np.linalg.svd(Wc[0].astype(np.float64), full_matrices=False)
rng = np.random.default_rng(7)
Q, _ = np.linalg.qr(rng.standard_normal((n, 32)))
addr = {"svd": Vt[:32].T.astype(np.float32), "rand": Q.astype(np.float32)}
edges = {1: 64, 2: 24, 3: 12}
# accumulators: per address type, per k: bin sums (B, n), counts (B,)
acc = {a: {k: (np.zeros((edges[k]**k, n)), np.zeros(edges[k]**k)) for k in edges} for a in addr}
S1 = np.zeros(n); S2 = np.zeros(n)
keep = 200_000; kF = np.empty((keep, n), dtype=np.float32); kA = {a: np.empty((keep, 32), dtype=np.float32) for a in addr}; kept = 0
B = 8192; done = 0; t0 = time.time()
def binid(t, k):  # t: (b, k) standard normal coords -> bin index with equiprobable edges
    nb = edges[k]; q = np.clip((norm_cdf(t) * nb).astype(int), 0, nb - 1)
    idx = np.zeros(len(t), dtype=int)
    for d in range(k): idx = idx * nb + q[:, d]
    return idx
from scipy.special import ndtr as norm_cdf
while done < T:
    b = min(B, T - done)
    x = rng.standard_normal((b, n), dtype=np.float32); h = x
    for l in range(L): h = np.maximum(h @ Wf[l], 0)
    F = h.astype(np.float64); S1 += F.sum(0); S2 += (F**2).sum(0)
    for a, P in addr.items():
        t = x @ P
        for k in edges:
            idx = binid(t[:, :k], k); sums, cnt = acc[a][k]
            np.add.at(sums, idx, F); np.add.at(cnt, idx, 1)
        if kept < keep: tk = min(b, keep - kept); kA[a][kept:kept + tk] = t[:tk, :32]
    if kept < keep: tk = min(b, keep - kept); kF[kept:kept + tk] = h[:tk]; kept += tk
    done += b
print(f"net {net}: {T} samples, {time.time()-t0:.0f}s", flush=True)
m = S1 / T; var = S2 / T - m**2; Vtot = var.sum()
print(f"total per-sample variance sum_i Var(F_i) = {Vtot:.4f} (avg {Vtot/n:.4f}); to match the chain at the 10% floor the residual must be <= {2.2e-8*6500/(Vtot/n):.4f} of it")
print("explained-variance fraction Var(E[F|T])/Var(F) by binning:")
for a in addr:
    for k in edges:
        sums, cnt = acc[a][k]; ok = cnt > 0; mb = sums[ok] / cnt[ok][:, None]
        expl = (cnt[ok][:, None] * (mb - m)**2).sum() / T
        print(f"  address {a:4s} k={k}: {expl/Vtot:.5f}  (bins {ok.sum()})")
# nearest-neighbour conditional variance for k = 8, 32 on the kept subsample
from scipy.spatial import cKDTree
kF = kF[:kept]
for a in addr:
    for k in (8, 32):
        A = kA[a][:kept, :k]; tree = cKDTree(A); d, j = tree.query(A, k=2); j = j[:, 1]
        cv = 0.5 * np.mean(np.sum((kF - kF[j])**2, axis=1))      # E Var(F|T) estimate (biased up by the NN distance)
        print(f"  address {a:4s} k={k} (NN pairs, mean NN distance {d[:,1].mean():.3f}): explained fraction <= {1 - cv/Vtot:.5f}")
