# Test of the "Ruelle bundle" (base space of 2^k sign chambers, Gaussian fibre per chamber) at its strongest:
# chamber probabilities and chamber-conditional means/covariances are taken EXACTLY (from samples), and we ask
# whether one Gaussian step per chamber, mixed with the chamber weights, has a smaller local defect than the
# single Gaussian step of the closure (same samples, same next-layer empirical mean).  The closure error is the
# sum of transported local defects (ledger), so this bounds what any base space can buy.
# Partitions of the source layer h_l (64 cells each):
#   W-frame : signs of the top-6 right singular vectors of the next weight matrix (the proposal)
#   data    : signs of the top-6 covariance eigenvectors of h_l (centred) -- best-case sign chambers
#   gain    : 64 quantile bins of |h_l|^2  (a one-dimensional base: the scale / gain coordinate)
#   random  : random labels (null: measures the sampling noise of the test)
import numpy as np, sys, time
from scipy.special import ndtr
from closure import phi
n, L, s, T = 256, 16, 0, int(float(sys.argv[1])) if len(sys.argv) > 1 else 1000000
SRC = [0, 2, 4, 6, 10, 14]                       # 0-based source layers h_l (predict layer l+2, 1-based)
Ws = np.load(f"W_n{n}_L{L}_s{s}.npy"); Wf = [W.T.astype(np.float32) for W in Ws]
hm = np.load(f"hmom_n{n}_L{L}_s{s}.npz")
P = {}
for l in SRC:
    U, S_, Vt = np.linalg.svd(Ws[l+1]); Pw = Vt[:6].astype(np.float32)
    C = hm["S2"][l] - np.outer(hm["S1"][l], hm["S1"][l]); ev, V = np.linalg.eigh(C); Pd = V[:, -6:].T.astype(np.float32)
    P[l] = (Pw, Pd, hm["S1"][l].astype(np.float32))
rng = np.random.default_rng(77); B = 20000
# pilot for gain quantiles
h = rng.standard_normal((B, n)).astype(np.float32); qs = {}
for l in range(L):
    h = np.maximum(h @ Wf[l], 0)
    if l in SRC: qs[l] = np.quantile((h*h).sum(1), np.linspace(0, 1, 65)[1:-1])
names = ["W-frame", "data", "gain", "random"]
cnt = {l: np.zeros((4, 64)) for l in SRC}; S1 = {l: np.zeros((4, 64, n)) for l in SRC}; S2 = {l: np.zeros((4, 64, n, n)) for l in SRC}
M1 = np.zeros((L, n)); bits = 2**np.arange(6)
t0 = time.time()
for b in range(T//B):
    h = rng.standard_normal((B, n)).astype(np.float32)
    for l in range(L):
        h = np.maximum(h @ Wf[l], 0); M1[l] += h.sum(0, dtype=np.float64)
        if l in SRC:
            Pw, Pd, mu = P[l]
            labs = [((h @ Pw.T) > 0) @ bits, (((h - mu) @ Pd.T) > 0) @ bits,
                    np.searchsorted(qs[l], (h*h).sum(1)), rng.integers(0, 64, B)]
            for p, lab in enumerate(labs):
                for c in range(64):
                    X = h[lab == c]
                    if len(X) == 0: continue
                    cnt[l][p, c] += len(X); S1[l][p, c] += X.sum(0, dtype=np.float64); S2[l][p, c] += X.T @ X
    if b % 10 == 9: print(f"  batch {b+1}/{T//B}  {time.time()-t0:.0f}s", flush=True)
M1 /= T
def mstep(mu, var):
    sig = np.sqrt(var); a = mu/sig; return sig*(a*ndtr(a) + phi(a))
print(f"n={n} L={L} s={s} T={T}: local defect of one Gaussian step (prediction of layer l+2 mean from moments of h_l)")
print(" layer | closure: rms   8*scale |  " + "  ".join(f"{nm:>8s}: MSE/closure 8*scale" for nm in names))
for l in SRC:
    W = Ws[l+1]; mstar = M1[l+1]; sc = lambda d: (mstar @ d)/(mstar @ mstar)
    N = cnt[l][0].sum(); m = S1[l][0].sum(0)/N; C = S2[l][0].sum(0)/N - np.outer(m, m)
    d0 = mstep(W @ m, np.einsum("ij,jk,ik->i", W, C, W)) - mstar
    row = f"  {l+2:3d}  | {np.sqrt(np.mean(d0**2)):.2e} {8*sc(d0):+.4f} |"
    for p in range(4):
        pred = np.zeros(n)
        for c in range(64):
            k = cnt[l][p, c]
            if k < 2: continue
            mc = S1[l][p, c]/k; Cc = S2[l][p, c]/k - np.outer(mc, mc)
            pred += (k/N)*mstep(W @ mc, np.einsum("ij,jk,ik->i", W, Cc, W))
        d = pred - mstar
        row += f"   {np.mean(d**2)/np.mean(d0**2):6.3f}     {8*sc(d):+.4f}   "
    print(row, flush=True)
