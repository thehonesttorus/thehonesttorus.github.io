# Post-hoc readout-weighted comparison for test S (note XXXVIII section 8): does the collective-symbol prediction R_pred
# capture the part of R that the next layer's mean map and covariance program read, compared with candidate D?
#   python symread.py NET LAYERS [MCPREFIX]      (needs symR_off{NET}_{l}.npz from symtest.py with SYMSAVE=1)
# Weights (the Hermite coefficients through which each slice of z_(l+1) enters, at the true reference of layer l+1):
#   diagonal: the mean map's kappa4 coefficient c4 = (alpha^2 - 1) phi(alpha) / (24 sigma^3);
#   (2,2) slice: its use in Cov(y_i, y_j), c(1,2)_i c(1,2)_j / 4 with c(1,2) = phi / sigma;
#   (3,1) slice: its use in Cov(y_i, y_j), c(1,3)_i c(1,1)_j / 6 with c(1,3) = -alpha phi / sigma^2, c(1,1) = Phi.
import sys, numpy as np
from math import erf, sqrt, pi
net = int(sys.argv[1]); LAY = [int(x) for x in sys.argv[2].split(",")]
pre = sys.argv[3] if len(sys.argv) > 3 else f"mc4_off{net}"
F = np.load(f"{pre}_full.npz")
Phi = np.vectorize(lambda x: 0.5 * (1 + erf(x / sqrt(2)))); phi = lambda x: np.exp(-x * x / 2) / sqrt(2 * pi)
en = lambda X: float(np.sum(X * X))
for l in LAY:
    z = np.load(f"symR_off{net}_{l}.npz"); ia, ib = z["ia"], z["ib"]
    mu = F["mu"][l + 1].astype(np.float64); sg = np.sqrt(F["var"][l + 1].astype(np.float64)); a = mu / sg
    c4 = (a * a - 1) * phi(a) / (24 * sg ** 3); c2 = phi(a) / sg; c3 = -a * phi(a) / sg ** 2; c1 = Phi(a)
    wts = [c4, c2[ia] * c2[ib] / 4, (c3[:, None] * c1[None, :]) / 6]
    R = [z[f"R{j}"] for j in range(3)]; D = [z[f"D{j}"] for j in range(3)]
    Ks = sorted({int(k[1:].split("_")[0]) for k in z.files if k.startswith("P")})
    out = [f"layer {l:2d} -> {l + 1:2d}  (readout-weighted explained share: diag | (2,2) | (3,1))"]
    def share(X):
        return [1 - en(w * (r - x)) / en(w * r) for w, r, x in zip(wts, R, X)]
    def fitw(X):     # coefficient fitted in the weighted metric
        return [float(np.vdot(w * r, w * x) / np.vdot(w * x, w * x)) for w, r, x in zip(wts, R, X)]
    bD = fitw(D); sD = share([b * x for b, x in zip(bD, D)])
    out.append("  D, coef fitted in this metric " + " ".join(f"{b:.2f}" for b in bD) + ": " + " ".join(f"{s:+.3f}" for s in sD))
    for K in Ks:
        P = [z[f"P{K}_{j}"] for j in range(3)]
        sP = share(P); bP = fitw(P)
        jt = []
        for w, r, x, d in zip(wts, R, P, D):
            Zm = np.stack([(w * x).ravel(), (w * d).ravel()], 1); cc, *_ = np.linalg.lstsq(Zm, (w * r).ravel(), rcond=None)
            jt.append(1 - en((w * r).ravel() - Zm @ cc) / en(w * r))
        out.append(f"  K={K:2d}: R_pred (coef 1) " + " ".join(f"{s:+.3f}" for s in sP) + "; coef in this metric "
                   + " ".join(f"{b:.2f}" for b in bP) + "; jointly with D " + " ".join(f"{s:+.3f}" for s in jt))
    print("\n".join(out), flush=True)
