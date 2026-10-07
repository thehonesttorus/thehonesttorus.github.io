# First-order gated transport of the legs (note XXXIV, source 1), tested against the measured omitted classes (note XXXVI).
#   python mclegs.py NET NSAMP SEED MCPREFIX
# Re-runs the Monte Carlo of `mcstats.py NET NSAMP SEED MCPREFIX post` (same seed, same batches, hence the same samples),
# with the gates Phi_l = Phi(mu_l / sigma_l) of the true pre-activation law taken from MCPREFIX_full.npz, and accumulates
# per layer, for U = X d (X = W_(l+1) o Phi_l, d = z_l - mu_l):
#     S_ai = E[d_a U_i^2],   P_ic = E[U_i^2 U_c],   u3_i = E[U_i^3].
# Their all-distinct parts (inclusion-exclusion with the true D21 and kappa3 diagonal of z_l, same samples) are the
# contractions of the true all-distinct third cumulant T of z_l that the theorem predicts, with no free coefficient:
#   kappa3, all-distinct y class transported (R3):  sum_(distinct) X_ia X_ib X_ic T_abc  (D3 slice), X_ia X_ib X_cd T_abd (D21)
#   kappa4, leg-fed (2,1,1) class transported (diag of R):  6 sum_a H_ia 2 m_a (1 - Phi_a) sum_(b != c, != a) X_ib X_ic T_abc
# The measured R3 and R (exact, same samples) are recomputed from MCPREFIX_full.npz as in yclasses.py, together with the
# Gaussian second-order candidates (one facet, two covariances) and the scale-mode shape D, and the fits are reported.
import sys, math, time, numpy as np
from math import lgamma, exp
net, N, seed, pre = int(sys.argv[1]), int(float(sys.argv[2])), int(sys.argv[3]), sys.argv[4]
Wcol = np.load(f"../official/W_off{net}.npy").astype(np.float32)
L, n, _ = Wcol.shape
F = np.load(f"{pre}_full.npz")
mu = F["mu"].astype(np.float64); var = F["var"].astype(np.float64); sd = np.sqrt(var); al = mu / sd
Phi = 0.5 * (1 + np.vectorize(math.erf)(al / math.sqrt(2))); phi = np.exp(-al * al / 2) / math.sqrt(2 * math.pi)
X32 = [(Wcol[l + 1] * Phi[l][None, :].astype(np.float32)) for l in range(L - 1)]
mu32 = mu.astype(np.float32)
B = 32768
rng = np.random.default_rng(seed)
S = np.zeros((L - 1, n, n)); P = np.zeros((L - 1, n, n)); u3 = np.zeros((L - 1, n)); cnt = 0
t0 = time.time(); done = 0; nb = 0
while done < N:
    h = 0 if done < N // 2 else 1
    b = min(B, N - done, (N // 2 - done) if h == 0 else N - done)
    Y = rng.standard_normal((n, b), dtype=np.float32)
    for l in range(L):
        Z = Wcol[l] @ Y
        if l < L - 1:
            D = Z - mu32[l][:, None]
            U = X32[l] @ D; U2 = U * U
            S[l] += D @ U2.T; P[l] += U2 @ U.T; u3[l] += (U2 * U).sum(1, dtype=np.float64)
        Y = np.maximum(Z, 0.0, out=Z)
    cnt += b; done += b; nb += 1
    if nb % 32 == 0:
        print(f"{done} samples, {time.time() - t0:.0f}s", flush=True)
S /= cnt; P /= cnt; u3 /= cnt
np.savez(f"{pre}_legs.npz", S=S.astype(np.float32), P=P.astype(np.float32), u3=u3)
print(f"accumulated {N} samples in {time.time() - t0:.0f}s", flush=True)

# ---------------- analysis ----------------
W64 = Wcol.astype(np.float64)
rng2 = np.random.default_rng(1)
ia = rng2.integers(0, n, 12000); ib = rng2.integers(0, n, 12000); k_ = ia != ib; ia, ib = ia[k_], ib[k_]


def offd(A):
    A = np.array(A, dtype=np.float64); np.fill_diagonal(A, 0.0); return A


def sym0(A):
    A = np.array(A, dtype=np.float64); A = 0.5 * (A + A.T); np.fill_diagonal(A, 0.0); return A


def T_pair(W, H, d=None, K=None, Bm=None):
    diag = np.zeros(n); k22 = np.zeros(ia.size); k31 = np.zeros((n, n)); W3 = W * H
    if d is not None:
        diag += (H * H) @ d; k22 += np.einsum("pa,pa->p", H[ia] * d, H[ib]); k31 += (W3 * d[None, :]) @ W.T
    if K is not None:
        HK = H @ K; Cab = W[ia] * W[ib]
        diag += 3 * np.einsum("ia,ia->i", HK, H)
        k22 += np.einsum("pa,pa->p", HK[ia], H[ib]) + 2 * np.einsum("pa,pa->p", Cab @ K, Cab)
        k31 += 3 * (HK * W) @ W.T
    if Bm is not None:
        G = Bm @ W.T; Xm = W * G.T; HX = H @ Xm.T
        diag += 4 * np.einsum("ia,ai->i", W3, G)
        k22 += 2 * (HX[ia, ib] + HX[ib, ia])
        k31 += W3 @ G + 3 * (H * G.T) @ W.T
    np.fill_diagonal(k31, 0.0)
    return diag, k22, k31


def T3_pair(W, H, k3v, Dm):
    G3 = Dm @ W.T
    d3 = (W * H) @ k3v + 3 * np.einsum("ia,ai->i", H, G3)
    d21 = (H * k3v[None, :]) @ W.T + H @ Dm @ W.T + 2 * (W * G3.T) @ W.T
    np.fill_diagonal(d21, 0.0)
    return d3, d21


def fitrep(R_, X_, fixed=True):
    """per slice: corr, explained with coefficient 1, fitted coefficient and its explained fraction"""
    out = []
    for r, x in zip(R_, X_):
        r, x = r.ravel(), x.ravel()
        c = float(np.corrcoef(r, x)[0, 1]); b = float(r @ x / (x @ x))
        e1 = 1 - float(np.linalg.norm(r - x) ** 2 / (r @ r)); eb = 1 - float(np.linalg.norm(r - b * x) ** 2 / (r @ r))
        out.append(f"corr {c:+.3f} coef1 {e1:+.3f} fit {b:+.3f}->{eb:+.3f}")
    return " | ".join(out)


def fit2(r, x1, x2):
    Z = np.stack([x1.ravel(), x2.ravel()], 1); c, *_ = np.linalg.lstsq(Z, r.ravel(), rcond=None)
    return c, 1 - float(np.linalg.norm(r.ravel() - Z @ c) ** 2 / (r.ravel() @ r.ravel()))


Er = math.sqrt(2.0 / n) * exp(lgamma((n + 1) / 2) - lgamma(n / 2)); GS = 2.0 / n; GM = (1 + 2.0 / n) - Er * Er * (n + 1.0) / n
print(f"net {net}: legs transport test on {pre} (same samples); slices kappa3 (D3 | D21 off-diag), kappa4 diag")
for l in range(L - 1):
    W = W64[l + 1]; H = W * W; X = W * Phi[l][None, :]; XX = X * X
    k3z = F["k3"][l].astype(np.float64); D21z = offd(F["D21"][l]); Cz = sym0(F["cov"][l])
    # all-distinct contractions of the true T of z_l
    tri_d3 = u3[l] - 3 * np.einsum("ia,ai->i", XX, D21z @ X.T) - (XX * X) @ k3z
    t1 = XX @ D21z @ X.T; t2 = 2 * (X * (D21z @ X.T).T) @ X.T; t3 = (XX * k3z[None, :]) @ X.T
    tri_d21 = offd(P[l] - t1 - t2 - t3)
    allD = S[l] - (XX @ D21z).T - 2 * X.T * (D21z @ X.T) - (XX * k3z[None, :]).T          # (a, i)
    # measured R3 and R (exact, same samples)
    my = F["mu_y"][l].astype(np.float64); vy = F["var_y"][l].astype(np.float64); k3y = F["k3_y"][l].astype(np.float64)
    Dy = offd(F["D21_y"][l]); Cy = F["cov_y"][l].astype(np.float64); Cy = 0.5 * (Cy + Cy.T); Cyo = offd(Cy)
    t3 = (F["k3"][l + 1].astype(np.float64), offd(F["D21"][l + 1]))
    R3 = tuple(a - b for a, b in zip(t3, T3_pair(W, H, k3y, Dy)))
    t4 = F["k4"][l + 1].astype(np.float64)
    R4d = t4 - T_pair(W, H, F["k4_y"][l].astype(np.float64), sym0(F["K22_y"][l]), offd(F["K31_y"][l]))[0]
    # candidates for R3: legs first order (coef 1), Gaussian second order (one facet, two covariances), diag only
    rho = phi[l] / sd[l]
    G = Cz @ X.T                                                     # G_ai = sum_b C_ab X_ib
    fac_d3 = 3 * np.einsum("ia,ia->i", W * rho[None, :], G.T ** 2 - XX @ (Cz * Cz))
    # candidates for R (kappa4 diag): leg-fed (2,1,1), and the scale-mode shape D (fitted)
    c21 = 2 * my * (1 - Phi[l])
    leg4 = 6 * np.einsum("ia,ai->i", H * c21[None, :], allD)
    mu1 = F["mu"][l + 1].astype(np.float64); S1 = F["cov"][l + 1].astype(np.float64); S1 = 0.5 * (S1 + S1.T); S1d = np.diag(S1).copy()
    k31 = F["k3"][l + 1].astype(np.float64); D1 = offd(F["D21"][l + 1])
    Dd = (GS * 3 * S1d ** 2 + GM * 4 * mu1 * k31) - T_pair(
        W, H, GS * 3 * vy * vy + GM * 4 * my * k3y,
        sym0(GS * (np.outer(vy, vy) + 2 * Cyo * Cyo) + GM * (2 * my[:, None] * Dy.T + 2 * my[None, :] * Dy)),
        offd(GS * 3 * vy[:, None] * Cyo + GM * (3 * my[:, None] * Dy + np.outer(k3y, my))))[0]
    rel = lambda a, b: float(np.linalg.norm(a) / np.linalg.norm(b))
    cDL, eDL = fit2(R4d, Dd, leg4)
    cLF, eLF = fit2(R3[0], tri_d3, fac_d3)
    eL1 = 1 - float(np.linalg.norm(R4d - leg4) ** 2 / (R4d @ R4d))
    resid = R4d - leg4; bD = float(resid @ Dd / (Dd @ Dd)); eDr = 1 - float(np.linalg.norm(resid - bD * Dd) ** 2 / (R4d @ R4d))
    print(f"layer {l:2d}->{l + 1:2d}: |R3|/|t3| {rel(R3[0], t3[0]):.3f} {rel(R3[1], t3[1]):.3f}; legs-1st-order |.|/|R3| "
          f"{rel(tri_d3, R3[0]):.3f} {rel(tri_d21, R3[1]):.3f}\n"
          f"    R3 vs legs (coef 1): {fitrep(R3, (tri_d3, tri_d21))}\n"
          f"    R3 D3 vs facet 2nd order: {fitrep((R3[0],), (fac_d3,))}; legs+facet fit coefs {cLF[0]:+.3f} {cLF[1]:+.3f} -> {eLF:.3f}\n"
          f"    R diag |R|/|t4| {rel(R4d, t4):.3f}; leg-fed |.|/|R| {rel(leg4, R4d):.3f}: {fitrep((R4d,), (leg4,))}; "
          f"D alone: {fitrep((R4d,), (Dd,))}\n"
          f"    R diag = a D + b leg-fed: a {cDL[0]:.3f} b {cDL[1]:+.3f} explained {eDL:.3f}; with leg-fed fixed at 1: "
          f"{eL1:+.3f}, then + D (coef {bD:.3f}): {eDr:+.3f}", flush=True)
