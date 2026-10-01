"""T19: share of the (2,1) slice D21_ab = kappa3(z_a, z_a, z_b) carried by its chaos-1 rank-one part sigma^2 (W^T t)_b
(row-constant), and the participation ratio of the slice before/after removing it. Own net n=256, depth 16, MC N."""
import sys, numpy as np
from common import bench
n, L, N = int(sys.argv[1]), 16, int(sys.argv[2]); chunk = 8192; s2 = 2.0 / n
W = bench.weights_from_seed(424242, n, L); layers = [5, 9, 13, 15]
def passes(fn):
    rng = np.random.default_rng(1)
    for c in range(N // chunk):
        h = rng.standard_normal((chunk, n)).astype(np.float32)
        for l in range(L):
            z = h @ W[l]; h = np.maximum(z, 0.0); fn(l, z, h)
sa = np.zeros((L, n)); sz = np.zeros((L, n))
def p1(l, z, h): sa[l] += h.sum(0, dtype=np.float64); sz[l] += z.sum(0, dtype=np.float64)
passes(p1); mu = sa / N; mz = sz / N
t = np.zeros((L, n)); D = {l: np.zeros((n, n)) for l in layers}
def p2(l, z, h):
    u = h.astype(np.float64) - mu[l]; t[l] += ((u ** 2).sum(1)[:, None] * u).sum(0)
    if l in layers:
        y = z.astype(np.float64) - mz[l]; D[l] += (y ** 2).T @ y
passes(p2); t /= N
def pr(M):
    sv = np.linalg.svd(M, compute_uv=False); return sv.sum() ** 2 / (sv ** 2).sum(), sv[0] ** 2 / (sv ** 2).sum()
for l in layers:
    Dl = D[l] / N; off = ~np.eye(n, dtype=bool)
    r1 = np.tile(s2 * (W[l].astype(np.float64).T @ t[l - 1]), (n, 1))
    Do = np.where(off, Dl, 0.0); Ro = np.where(off, r1, 0.0); Res = Do - Ro
    e = 1 - (Res ** 2).sum() / (Do ** 2).sum()
    print(f"z layer {l+1}: off-diag energy explained by rank-one trace part {e:.3f}; PR(slice) {pr(Do)[0]:.1f} (top-1 share {pr(Do)[1]:.3f}) "
          f"-> PR(after) {pr(Res)[0]:.1f} (top-1 {pr(Res)[1]:.3f})", flush=True)
    U, sv, Vt = np.linalg.svd(Res); u1, v1 = U[:, 0], Vt[0]
    Wl = W[l].astype(np.float64); m = mz[l]; var = np.diag((D[l] * 0)) if False else None
    cands_r = dict(m=m, absm=np.abs(m), one=np.ones(n), Wt_t=Wl.T @ t[l - 1], Wt_mu=Wl.T @ mu[l - 1])
    cands_c = dict(m=m, one=np.ones(n), Wt_t=Wl.T @ t[l - 1], Wt_mu=Wl.T @ mu[l - 1])
    cr = lambda a, b: abs(np.corrcoef(a, b)[0, 1]) if np.std(b) > 0 else abs(a.mean()) / np.sqrt((a ** 2).mean())
    print("   left u1 vs", {k: round(cr(u1, v), 3) for k, v in cands_r.items()}, " right v1 vs", {k: round(cr(v1, v), 3) for k, v in cands_c.items()})
