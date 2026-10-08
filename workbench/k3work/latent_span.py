# Column-space test of the latent (collective) model for the true slices, in consistent bases.
#   python latent_span.py NET   (mc2_off{NET}_full.npz)
# A K-mode latent model z - mu = U t + (independent rest) makes every off-diagonal slice linear in U_j on its single
# index: K31_ij = T4[U_i, U_i, U_i, U_j], D21_ij = T3[U_i, U_i, U_j]. So the plain slice's column space must lie in
# span(U) (U = top-K eigenvectors of C), and the read-weighted slice X diag(Phi) in span(diag(Phi) U). Reported for the
# truth and for its gain residual X - g S (g fitted plainly), at K = 16, 64, 256; and, for comparison, the share in a
# random K-dimensional subspace (K / n).
import sys, numpy as np
from math import sqrt, pi, erf
net = int(sys.argv[1])
T0 = np.load(f"mc2_off{net}_full.npz"); T = {k: T0[k] for k in T0.files}
phi = lambda a: np.exp(-0.5 * a * a) / sqrt(2 * pi)
Phi = np.vectorize(lambda a: 0.5 * (1 + erf(a / sqrt(2))))
f64 = lambda a: np.asarray(a, np.float64)
def offd(M):
    M = M.copy(); np.fill_diagonal(M, 0.0); return M
def colshare(X, Q):
    return float(np.sum((X @ Q) ** 2)) / float(np.sum(X * X))
for l in range(2, 16):
    mu, var = f64(T["mu"][l]), f64(T["var"][l]); C = f64(T["cov"][l]); s = np.sqrt(var); al = mu / s
    c1 = Phi(al)
    ev, U = np.linalg.eigh(C); U = U[:, ::-1]
    D21, K31 = offd(f64(T["D21"][l])), offd(f64(T["K31"][l]))
    S21 = offd(0.5 * (2 * mu[:, None] * C + var[:, None] * mu[None, :])); S31 = offd(3 * var[:, None] * C)
    line = [f"layer {l:2d}"]
    for nm, X, S in (("D21", D21, S21), ("K31", K31, S31)):
        g = float(np.sum(X * S)) / float(np.sum(S * S)); R = X - g * S
        parts = []
        for K in (16, 64, 256):
            Q = U[:, :K]; Qg, _ = np.linalg.qr(c1[:, None] * Q)
            parts.append(f"K{K}: plain {colshare(X, Q):.3f}/res {colshare(R, Q):.3f}; read {colshare(X * c1[None, :], Qg):.3f}/res {colshare(R * c1[None, :], Qg):.3f}")
        line.append(f"  {nm}: " + "  ".join(parts))
    print("\n".join(line), flush=True)
print("random K-dim subspace share: K/n =", ", ".join(f"{K / 1024:.3f}" for K in (16, 64, 256)))
