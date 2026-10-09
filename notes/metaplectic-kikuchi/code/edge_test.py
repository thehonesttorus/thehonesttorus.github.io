import sys, numpy as np
sys.argv = ["x", ".", ".", "0", "1"]
src = open("/home/user/thehonesttorus.github.io/workbench/k3work/var_ladder.py").read()
src = src.split("cd = np.load")[0]          # definitions only
exec(src)
rng = np.random.default_rng(0)
n = 4; N = 4_000_000
A = rng.standard_normal((n, n)) * 0.5 + np.eye(n)
mu = np.array([0.3, -0.5, 1.0, 0.0])
for eps in (0.08, 0.04):
    Z = rng.standard_normal((N, n))
    # non-Gaussian: quadratic perturbation with cross terms
    Q = rng.standard_normal((n, n)) * 0.5
    H = mu + Z @ A.T + eps * ((Z ** 2 - 1) @ Q.T) + eps * 0.7 * (Z[:, [0]] * Z[:, [1]]) * np.array([1, -1, 0.5, 0.3])
    X = np.maximum(H, 0)
    covX = np.cov(X.T)
    m = H.mean(0); Hc = H - m
    C = Hc.T @ Hc / N
    k3 = (Hc ** 3).mean(0); k4 = (Hc ** 4).mean(0) - 3 * np.diag(C) ** 2
    D21 = (Hc ** 2).T @ Hc / N
    K22 = (Hc ** 2).T @ (Hc ** 2) / N - np.outer(np.diag(C), np.diag(C)) - 2 * C ** 2
    K31 = (Hc ** 3).T @ Hc / N - 3 * np.diag(C)[:, None] * C
    KGm = KG(m, C)
    D3, D4 = edgeworth(m, C, k3, k4, D21, K22, K31)
    err0 = covX - KGm; err3 = covX - KGm - D3; err4 = covX - KGm - D3 - D4
    print(f"eps {eps}: |cov - KG| {np.abs(err0).max():.2e}  after E3 {np.abs(err3).max():.2e}  after E3+E4 {np.abs(err4).max():.2e}  (MC se ~{np.sqrt(np.diag(covX)).max()**2/np.sqrt(N)*2:.1e})")
    print("   diag:", np.round(np.diag(err0), 5), np.round(np.diag(err4), 5))
    print("   off :", np.round(err0[0, 1:], 5), np.round(err4[0, 1:], 5))
