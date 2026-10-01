"""layer-by-layer check of an estimator's pre-activation statistics against Monte Carlo (width-64 MLP seed 0)."""
import numpy as np, stageq, sys, importlib, os
mod = importlib.import_module(sys.argv[1]); N = float(sys.argv[2]); L = int(sys.argv[3]) if len(sys.argv) > 3 else 6
cfg = eval(sys.argv[4]) if len(sys.argv) > 4 else {}
W, mean, noise, _ = stageq.load(64, 0); n = 64
est, dg = mod.estimate(W, cfg, return_diag=True)
cache = f"/root/sg/mcstats_w64s0_{int(N)}.npz"
if os.path.exists(cache):
    d = np.load(cache); S, C, C21, done = d['S'], d['C'], d['C21'], int(d['done'])
else:
    rng = np.random.default_rng(5)
    S = np.zeros((8, 4, n)); C = np.zeros((8, n, n)); C21 = np.zeros((8, n, n)); done = 0
    mu_true = None
    while done < N:
        x = rng.standard_normal((1 << 17, n)); a = x
        for l in range(8):
            z = a @ W[l]; a = np.maximum(z, 0)
            zc = z - mean[l - 1] @ W[l] if l > 0 else z     # centre with near-exact means
            for k in range(4): S[l, k] += (zc ** (k + 1)).sum(0)
            C[l] += zc.T @ zc; C21[l] += (zc ** 2).T @ zc
        done += 1 << 17
    np.savez(cache, S=S, C=C, C21=C21, done=done)
off = ~np.eye(n, dtype=bool)
for l in range(1, L):
    m = S[l] / done; mu = m[0]; v = m[1] - mu ** 2
    c3 = m[2] - 3 * mu * m[1] + 2 * mu ** 3
    k4 = m[3] - 4 * mu * m[2] + 6 * mu ** 2 * m[1] - 3 * mu ** 4 - 3 * v ** 2
    Cz = C[l] / done - np.outer(mu, mu)
    E2 = C[l] / done; E21 = C21[l] / done
    D21 = E21 - np.outer(m[1], mu) - 2 * E2 * mu[None, :] + 2 * np.outer(mu ** 2, mu)  # kappa(z_a,z_a,z_c), rows a
    d = dg[l - 1]
    r = lambda x: np.sqrt(np.mean(x ** 2))
    s = f"z_{l+1}: var rel {r(d['vz']/v-1):.1e} offcov {r((d['Cz']-Cz)[off]):.1e}/{r(Cz[off]):.1e} k3 {r(d['k3']-c3):.1e}/{r(c3):.1e} k4 {r(d['k4']-k4):.1e}/{r(k4):.1e}"
    if 'D' in d: s += f"  D21off {r((d['D']-D21)[off]):.1e}/{r(D21[off]):.1e}"
    s += f"  mean(a) mse {np.mean((est[l]-mean[l])**2):.1e}"
    print(s)
