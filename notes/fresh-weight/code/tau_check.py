import sys, numpy as np
mcd, cdd, net = sys.argv[1], sys.argv[2], int(sys.argv[3])
c2 = np.load(f"{cdd}/chaindump2_{net}.npz"); F = np.load(f"{mcd}/mc2_off{net}_full.npz")
cov, mu = F["cov"], F["mu"]
print(f"net {net}: layer  tau_true = lam1/(u.mu)^2   measured dlam1/lam1 (u^T E u / lam1)   -tau/2   ratio")
for s in range(2, 15):
    Ct = 0.5 * (cov[s].astype(np.float64) + cov[s].astype(np.float64).T)
    Cc = 0.5 * (c2[f"C_off_{s}"].astype(np.float64) + c2[f"C_off_{s}"].astype(np.float64).T); np.fill_diagonal(Cc, c2[f"var_{s}"])
    lam, V = np.linalg.eigh(Ct); u = V[:, -1]; l1 = lam[-1]; m = mu[s].astype(np.float64)
    tau = l1 / (u @ m) ** 2
    d = float(u @ (Cc - Ct) @ u) / l1
    print(f"   {s:2d}   {tau:.4e}   {d:+.3e}   {-tau / 2:+.3e}   {d / (-tau / 2):.2f}")
