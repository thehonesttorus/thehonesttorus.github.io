# Is the collective-bulk cross part of the fourth-cumulant slices a shared matrix gain?
#   python latent_mgain.py NET "CHUNK_GLOB"   (mclatent.py chunks; reuses latent_an.py's merge and cumulant formulas)
# Model: z - mu = U t + r, r | t ~ N(0, diag(s^2) (1 + h(t))) with one modulation h(t) = beta . t + t^T G t - E[...]
# shared by all units (the scalar gain is G proportional to the inverse latent covariance). Then, for i != j,
#   K31_ij = T4[U_i, U_i, U_i, U_j] + 6 s_i^2 U_i^T M U_j + 3 s_i^2 T3[U_i, beta, U_j],   M = Lam G Lam,
#   k4_i  = T4[U_i^4] + 12 s_i^2 U_i^T M U_i + 3 s_i^4 Var h + (beta terms),
# so the cross part X = R_d - R_c - R_r of the non-gain (3,1) slice (latent_an.py) should be diag(s^2) U M' U^T with
# one K x K matrix M' (= 6 M plus the beta term), and the diagonal's cross part 2 s_i^2 U_i^T M' U_i with the same M'.
# Fitted by least squares (M' = (U^T S^2 U)^-1 U^T S X U), against the upper bound of any row-wise collective matrix
# (X U U^T); on half 0 and scored on half 1 as well, so that K^2 parameters are not credited with noise.
import sys, glob, re, numpy as np
from math import sqrt, pi, erf
sys.argv = sys.argv[:3]
exec(open("latent_an.py").read().split("chunks = {")[0])          # imports, load(), finish()
net, pat = int(sys.argv[1]), sys.argv[2]
files = sorted(glob.glob(pat)); chunks = {int(re.search(r"_(\d+)\.npz$", f).group(1)): f for f in files}
full = load(files); h0 = load([f for i, f in chunks.items() if i % 2 == 0]); h1 = load([f for i, f in chunks.items() if i % 2 == 1])
K = full["K"]; T = np.load(f"mc2_off{net}_full.npz")
en = lambda X: float(np.sum(X * X))
print(f"net {net}: {int(full['c'])} samples, K = {K}")


def parts(a, l):
    F = {v: finish(a, l, v) for v in ("d", "c", "r")}
    S = {v: 3 * F[v]["var"][:, None] * F[v]["cov"] for v in F}
    for v in S:
        np.fill_diagonal(S[v], 0.0)
    g = float(np.sum(F["d"]["K31"] * S["d"])) / en(S["d"])
    Rd, Rc, Rr = (F[v]["K31"] - g * S[v] for v in ("d", "c", "r"))
    X = Rd - Rc - Rr
    g4 = float(np.sum(F["d"]["k4"] * 3 * F["d"]["var"] ** 2)) / en(3 * F["d"]["var"] ** 2)
    R4 = {v: F[v]["k4"] - g4 * 3 * F[v]["var"] ** 2 for v in F}
    X4 = R4["d"] - R4["c"] - R4["r"]
    return F, Rd, Rc, X, X4, R4["d"]


for l in [int(x) for x in full["layers"]]:
    ev, U = np.linalg.eigh(np.asarray(T["cov"][l], np.float64)); U = U[:, ::-1][:, :K]
    F, Rd, Rc, X, X4, R4d = parts(full, l)
    _, _, _, X0, X40, _ = parts(h0, l); _, Rd1, Rc1, X1, X41, R4d1 = parts(h1, l)
    s2 = F["r"]["var"]
    def fitM(Xf, Uk):
        A = Uk.T @ (s2[:, None] ** 2 * Uk)
        return np.linalg.solve(A, Uk.T @ (s2[:, None] * Xf) @ Uk)
    out = [f"layer {l:2d}"]
    for k in (8, 16, 32, K):
        Uk = U[:, :k]; Mf = fitM(X, Uk); M0 = fitM(X0, Uk)
        pred = s2[:, None] * (Uk @ Mf @ Uk.T); np.fill_diagonal(pred, 0.0)
        p0 = s2[:, None] * (Uk @ M0 @ Uk.T); np.fill_diagonal(p0, 0.0)
        ub = X @ Uk @ Uk.T
        d4 = 2 * s2 * np.einsum("ip,pq,iq->i", Uk, 0.5 * (Mf + Mf.T), Uk)   # the diagonal's prediction from the same M'
        c4 = float(np.dot(d4, X4)) / max(np.linalg.norm(d4) * np.linalg.norm(X4), 1e-300)
        out.append(f"  K={k:2d}: cross part of the (3,1) slice explained by shared M' {1 - en(X - pred) / en(X):+.3f} (fit on half 0, scored on"
                   f" half 1: {1 - en(X1 - p0) / en(X1):+.3f}); any row-wise collective matrix {1 - en(X - ub) / en(X):+.3f};"
                   f" diagonal cross part: corr with 2 s^2 U^T M' U {c4:+.3f}, amplitude {float(np.dot(d4, X4)) / en(d4):+.2f}")
    Uk = U[:, :K]; Mf = fitM(X, Uk); pred = s2[:, None] * (Uk @ Mf @ Uk.T); np.fill_diagonal(pred, 0.0)
    tot = Rc + pred
    out.append(f"  non-gain (3,1) slice: collective T4 alone {1 - en(Rd - Rc) / en(Rd):+.3f}; + shared matrix gain (K={K}) "
               f"{1 - en(Rd - tot) / en(Rd):+.3f}; cross share of it {en(X) / en(Rd):.3f}; M' symmetric part share "
               f"{en(0.5 * (Mf + Mf.T)) / en(Mf):.3f}, top eigenvalues of sym M' " + " ".join(f"{x:+.2e}" for x in np.linalg.eigvalsh(0.5 * (Mf + Mf.T))[::-1][:4]))
    print("\n".join(out), flush=True)
