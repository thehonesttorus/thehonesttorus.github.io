import sys, numpy as np
sys.path.insert(0, '../../bench'); import bench
S = bench.load_set(sys.argv[1]); W = bench.weights(S, int(sys.argv[2])).astype(np.float64); N = int(float(sys.argv[3]))
L, n, _ = W.shape
acc = [dict() for _ in range(L)]; rng = np.random.default_rng(3); done = 0; B = 100000
sphere = len(sys.argv) > 4
while done < N:
    a = rng.standard_normal((B, n))
    if sphere: a *= np.sqrt(n) / np.linalg.norm(a, axis=1, keepdims=True)
    for l in range(L):
        z = a @ W[l]; a = np.maximum(z, 0); d = acc[l]
        for k, v in (('1', z.sum(0)), ('2', z.T @ z), ('22', (z * z).T @ (z * z)), ('21', (z * z).T @ z), ('3', (z**3).sum(0)), ('4', (z**4).sum(0))):
            d[k] = d.get(k, 0) + v
    done += B
for l in [1, 3, 6, 10, 15]:
    d = {k: v / N for k, v in acc[l].items()}
    m = d['1']; zc = None
    # central moments via raw: use cumulant formula for kappa(a,a,b,b)
    E1 = m; E2 = d['2']; E21 = d['21']; E22 = d['22']
    C = E2 - np.outer(m, m); Ea2 = np.diag(E2)
    # center: compute kappa22 = E[(a-ma)^2 (b-mb)^2] - C_aa C_bb - 2 C_ab^2
    ma, mb = m[:, None], m[None, :]
    E12 = E21.T
    c22 = (E22 - 2 * mb * E21 - 2 * ma * E12 + mb**2 * Ea2[:, None] + ma**2 * Ea2[None, :] + 4 * ma * mb * E2
           - 3 * ma**2 * mb**2 - 2 * ma * mb**2 * 0 )  # placeholder, recompute exactly below
    # exact: E[(a-ma)^2 (b-mb)^2] = E22 - 2mb E21 - 2ma E12 + mb^2 Ea2 + ma^2 Eb2 + 4 ma mb E11 - 3 ma^2 mb^2
    c22 = E22 - 2 * mb * E21 - 2 * ma * E12 + mb**2 * Ea2[:, None] + ma**2 * Ea2[None, :] + 4 * ma * mb * E2 - 3 * ma**2 * mb**2
    Q = c22 - np.outer(np.diag(C), np.diag(C)) - 2 * C**2
    off = ~np.eye(n, dtype=bool)
    Qo = Q.copy(); np.fill_diagonal(Qo, 0)
    ev = np.linalg.eigvalsh(Q)[::-1]
    v = np.diag(C)
    g = np.outer(v, v); coef = (Q[off] * g[off]).sum() / (g[off] ** 2).sum()
    resid = Q[off] - coef * g[off]
    print(f"l={l}: rms Q off {np.sqrt(np.mean(Q[off]**2)):.4f} mean {Q[off].mean():.4f} | top eig {ev[:4].round(3)} sum|ev| {np.abs(ev).sum():.2f} | "
          f"Q ~ eps v v^T: eps {coef:.4f} (2/n={2/n:.4f}), resid frac {np.sqrt(np.mean(resid**2))/np.sqrt(np.mean(Q[off]**2)):.3f}", flush=True)
