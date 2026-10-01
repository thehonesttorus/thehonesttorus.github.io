"""Is the rank-one spike of Q = kappa(z_a,z_a,z_b,z_b) carried by the top principal coordinate g of z?"""
import sys, numpy as np
sys.path.insert(0, '../../bench'); import bench
S = bench.load_set(sys.argv[1]); W = bench.weights(S, int(sys.argv[2])).astype(np.float64); N = int(float(sys.argv[3]))
L, n, _ = W.shape
rng = np.random.default_rng(5)
x = rng.standard_normal((N, n)); a = x
def Qmat(Z):
    Zc = Z - Z.mean(0); C = Zc.T @ Zc / len(Z); Z2 = Zc * Zc
    return Z2.T @ Z2 / len(Z) - np.outer(np.diag(C), np.diag(C)) - 2 * C ** 2, C
off = ~np.eye(n, dtype=bool)
for l in range(L):
    z = a @ W[l]; a = np.maximum(z, 0)
    if l in (1, 3, 6, 10, 15):
        Q, C = Qmat(z)
        ev, U = np.linalg.eigh(C); u = U[:, -1]
        g = (z - z.mean(0)) @ u; g = g / g.std()
        ell = ((z - z.mean(0)) * g[:, None]).mean(0)
        k3g = (g ** 3).mean(); k4g = (g ** 4).mean() - 3
        # residual after removing the linear dependence on g
        r = z - z.mean(0) - np.outer(g, ell)
        Qr, Cr = Qmat(r)
        # heteroscedastic slopes: E[r^2 | g] ~ v + s g + t (g^2-1)
        r2 = r * r
        s = (r2 * g[:, None]).mean(0); t = (r2 * (g[:, None] ** 2 - 1)).mean(0) / 2
        Qhet = np.outer(s, s) + 2 * np.outer(t, t)
        rms = lambda M: np.sqrt(np.mean(M[off] ** 2))
        resid = Qr - Qhet
        print(f"l={l}: top-eig frac {ev[-1]/ev.sum():.3f} k3(g) {k3g:.3f} k4(g) {k4g:.3f} | rms Q {rms(Q):.4f} mean {Q[off].mean():.4f} | "
              f"residual Q rms {rms(Qr):.4f} mean {Qr[off].mean():.4f} | hetero-slope model rms {rms(Qhet):.4f}, Qr - model rms {rms(resid):.4f} mean {resid[off].mean():.4f}", flush=True)
