# Rows of the true (3,1) slice against rows of D21 and C (note XXXIX). The per-row least-squares coefficients c_a of
# K31[a, :] on D21[a, :] (and d_a on C[a, :]) explain much of K31 while a single coefficient explains nothing; this asks
# whether c_a is a function of the unit's own carried state, and what a two-regressor row model (C and D21) explains.
#   python k31rows.py NET MCFULL MCH0 MCH1
import sys, numpy as np

net, ff, f0, f1 = int(sys.argv[1]), sys.argv[2], sys.argv[3], sys.argv[4]
F, H0, H1 = np.load(ff), np.load(f0), np.load(f1)
n = F["mu"].shape[1]


def zd(x):
    x = np.array(x, dtype=np.float64); np.fill_diagonal(x, 0.0); return x


print(f"net {net}: two-regressor row model K31[a,:] ~ d_a C[a,:] + c_a D21[a,:] (explained share of noise-corrected energy); "
      "correlation of c_a (D21 coefficient, from the joint fit) with per-unit features; and the share explained when c_a, d_a "
      "are replaced by their best fit on the features")
for l in range(2, 16):
    K = zd(F["K31"][l].T); dn = zd(H0["K31"][l].T) - zd(H1["K31"][l].T)
    E = np.sum(K ** 2) - np.sum(dn ** 2) / 4
    var = F["var"][l].astype(np.float64); mu = F["mu"][l].astype(np.float64); sd = np.sqrt(var); al = mu / sd
    k3 = F["k3"][l].astype(np.float64); k4 = F["k4"][l].astype(np.float64)
    C = zd(F["cov"][l]); D = zd(F["D21"][l])
    # joint per-row fit
    cc = np.sum(C * C, 1); dd = np.sum(D * D, 1); cd = np.sum(C * D, 1); kc = np.sum(K * C, 1); kd = np.sum(K * D, 1)
    det = np.maximum(cc * dd - cd * cd, 1e-300)
    dco = (kc * dd - kd * cd) / det; cco = (kd * cc - kc * cd) / det
    expl2 = np.sum(dco * kc + cco * kd) / E
    expl_c = np.sum(kc ** 2 / np.maximum(cc, 1e-300)) / E
    feats = {"alpha": al, "1/alpha": 1 / np.where(np.abs(al) > 1e-3, al, 1e-3), "var/mu": var / np.where(np.abs(mu) > 1e-9, mu, 1e-9),
             "k3/var^1.5": k3 / var ** 1.5, "k4/k3": k4 / np.where(np.abs(k3) > 1e-12, k3, 1e-12), "sd": sd, "mu": mu}
    # weight units by how much they matter in the fit (|D21 row|^2 |c|...): plain Pearson over units with weight dd
    wts = dd / dd.sum()
    def wcorr(x, y):
        x = x - np.sum(wts * x); y = y - np.sum(wts * y)
        return np.sum(wts * x * y) / np.sqrt(np.sum(wts * x * x) * np.sum(wts * y * y) + 1e-300)
    cors = " ".join(f"{k} {wcorr(v, cco * sd):+.2f}" for k, v in feats.items())
    # feature model: c_a sd_a ~ polynomial in alpha (degree 3), d_a ~ polynomial in alpha; refit globally
    A = np.stack([al ** p for p in range(4)], 1)
    X = np.concatenate([A[:, :, None] * (C / 1.0)[:, None, :], A[:, :, None] * (D / sd[:, None])[:, None, :]], axis=1)  # n x 8 x n
    G = np.einsum("apn,aqn->pq", X, X); g = np.einsum("apn,an->p", X, K)
    b = np.linalg.solve(G + 1e-12 * np.trace(G) * np.eye(len(g)), g)
    explf = (2 * b @ g - b @ G @ b) / E
    print(f"{l:2d} | C rows {expl_c:.3f}, C + D21 rows {expl2:.3f}, alpha-polynomial model {explf:.3f} | c_a sd_a vs {cors}", flush=True)
