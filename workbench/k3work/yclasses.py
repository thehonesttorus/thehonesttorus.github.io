# The one-step kappa4 error of the chain split by the support classes of the post-activation tensor (note XXXVI).
#   python yclasses.py NET MCPREFIX [CHAINPREFIX]
# MCPREFIX_{full,h0,h1}.npz from `mcstats.py NET N SEED MCPREFIX post`: slices of z_l and of y_l = relu(z_l) on the same
# samples. The empirical cumulant tensor of z_(l+1) = W y_l is exactly W# of the empirical cumulant tensor of y_l, so on any
# one sample
#     t = T_pair(y) + R,   T_pair(y) = W#[D(k4_y) + S22(K22_y) + S31(K31_y)]   (exact, (2,2) cross term included),
# and R is exactly the transport of the y (2,1,1) and (1,1,1,1) classes (the omitted classes of notes XXXIV and XXXV).
# Reported per transition l -> l+1, slices diag | (2,2) on random pairs | (3,1) off-diagonal:
#   - |R| / |t|, the split-half correlation of R and its noise share;
#   - the closed-form candidates for R that the theory names, each fitted with one scalar per layer:
#       A  the chain's lam core (diag 2 so2, (2,2) (so2_i + so2_j)/3, (3,1) C_off(z')), so2 = Var z' - H Var y;
#       B  the scale mode: g Sym(C_y x C_y) transported, minus its own pair-supported part (y-centred covariance);
#       C  the Gaussian-reference second-order triple class (products of two z covariances, gated), diag only, whose
#          coefficient is fixed at 1 by the theory;
#       D  the input-radius mode, parameter-free: for a bias-free ReLU net with x ~ N(0, I_n), z_l = r zeta_l exactly
#          (r = |x|/sqrt(n) independent of the direction), so by total cumulance kappa4 carries
#          Var(r^2) Sym3(Sigma, Sigma) + Cov(r, r^3) Sym4(mu, kappa3) at every layer (Var(r^2) = 2/n, Cov(r, r^3) ~ 3/(2n));
#          W# maps it to the same form in z' (mu, Sigma, kappa3 of z'), so its omitted-class part is the transported
#          full form minus T_pair of its own pair-supported y part. Reported with coefficient 1 and fitted. (Theory note:
#          in the Gaussian-input process this part should cancel against the sub-Gaussian kappa4 that the direction
#          zeta inherits from the sphere, -2/(n+2) Sym3(I, I) at the input, so its coefficient should be ~0. A fitted
#          coefficient measures a raw-vector scale mode of any origin);
#   - the gain-mode prediction for the scale coefficients, parameter-free from pair-supported y statistics:
#       g_raw = Var|y|^2 / (E|y|^2)^2 and g_cen = Var|y - m|^2 / (E|y - m|^2)^2, each also minus the radial 2/n;
#       g_conn = connected part of Var|y|^2 over (E|y|^2)^2 = [sum_ab kappa(y_a,y_a,y_b,y_b) + 4 sum_ab m_a kappa(y_a,y_b,y_b)] / (E|y|^2)^2,
#       the excess scale variance of y_l (the trace sector). Prediction: R ~ (g_conn / Var(r^2)) D, i.e. candidate D's shape
#       with its coefficient fixed by the trace sector (any scale fluctuation of y_l is transported downstream exactly,
#       by homogeneity, with the radial-tangent shape; the input radius itself contributes nothing net);
#   - with CHAINPREFIX (chain_off{NET}_o.npz, the one-step dump whose oracle inputs were CHAINPREFIX_h0): the pair
#     program's closure error on the y slices it computes (K4v, K22) and of note XXXIV's (3,1)_y closed form, each
#     transported, and the exact split of the chain's slice error
#         own - t = [c6x - T_pair(y)] + [lam core - R] + [own - c6x - lam core]
#     (closure | omitted classes net of the lam core | representation), c6x = exact transport of the chain's y slices
#     with B^y from the closed form.
import sys, os, math, numpy as np
from math import lgamma, exp
net, pre = int(sys.argv[1]), sys.argv[2]
cpre = sys.argv[3] if len(sys.argv) > 3 else None
Wcol = np.load(f"../official/W_off{net}.npy").astype(np.float64)
L, n, _ = Wcol.shape
M = {h: np.load(f"{pre}_{h}.npz") for h in ("full", "h0", "h1")}
ch = np.load(f"chain_off{net}_o.npz") if cpre and os.path.exists(f"chain_off{net}_o.npz") else None
I0 = np.load(f"{cpre}_h0.npz") if ch is not None else None
Phi = lambda x: 0.5 * (1 + np.vectorize(math.erf)(x / math.sqrt(2))); phi = lambda x: np.exp(-x * x / 2) / math.sqrt(2 * math.pi)
rng = np.random.default_rng(1)
n_in = n   # input dimension (all widths are n)
Er = math.sqrt(2.0 / n_in) * exp(lgamma((n_in + 1) / 2) - lgamma(n_in / 2))
GS = 2.0 / n_in                                            # Var(r^2), r^2 = chi^2_n / n
GM = (1 + 2.0 / n_in) - Er * Er * (n_in + 1.0) / n_in      # Cov(r, r^3) = E r^4 - E r E r^3
ia = rng.integers(0, n, 12000); ib = rng.integers(0, n, 12000); k_ = ia != ib; ia, ib = ia[k_], ib[k_]


def offd(X):
    X = np.array(X, dtype=np.float64); np.fill_diagonal(X, 0.0); return X


def sym0(X):
    X = np.array(X, dtype=np.float64); X = 0.5 * (X + X.T); np.fill_diagonal(X, 0.0); return X


def T_pair(W, H, d=None, K=None, B=None):
    """exact slices of W#[D(d) + S22(K) + S31(B)]: diag, (2,2) on the pairs (ia, ib), (3,1)_ij = kappa(i,i,i,j) (zero diag).
    K symmetric zero-diagonal, B_ab = kappa(y_a, y_a, y_a, y_b) zero-diagonal."""
    diag = np.zeros(n); k22 = np.zeros(ia.size); k31 = np.zeros((n, n)); W3 = W * H
    if d is not None:
        diag += (H * H) @ d; k22 += np.einsum("pa,pa->p", H[ia] * d, H[ib]); k31 += (W3 * d[None, :]) @ W.T
    if K is not None:
        HK = H @ K; Cab = W[ia] * W[ib]
        diag += 3 * np.einsum("ia,ia->i", HK, H)
        k22 += np.einsum("pa,pa->p", HK[ia], H[ib]) + 2 * np.einsum("pa,pa->p", Cab @ K, Cab)
        k31 += 3 * (HK * W) @ W.T
    if B is not None:
        G = B @ W.T; X = W * G.T; HX = H @ X.T
        diag += 4 * np.einsum("ia,ai->i", W3, G)
        k22 += 2 * (HX[ia, ib] + HX[ib, ia])
        k31 += W3 @ G + 3 * (H * G.T) @ W.T
    np.fill_diagonal(k31, 0.0)
    return diag, k22, k31


def truth(S, l):
    return S["k4"][l + 1].astype(np.float64), S["K22"][l + 1].astype(np.float64)[ia, ib], offd(S["K31"][l + 1])


def ytens(S, l):
    return S["k4_y"][l].astype(np.float64), sym0(S["K22_y"][l]), offd(S["K31_y"][l])


def add(*Xs):
    return tuple(sum(z) for z in zip(*Xs))


def sub(X, Y):
    return tuple(a - b for a, b in zip(X, Y))


def nm(X):
    return [float(np.linalg.norm(x)) for x in X]


def gate(mu, var):
    s = np.sqrt(var); a = mu / s; Ph = Phi(a); ph = phi(a)
    m = mu * Ph + s * ph; v = (mu * mu + var) * Ph + mu * s * ph - m * m
    return dict(s=s, a=a, Ph=Ph, ph=ph, m=m, v=v, dl=ph / s, dl1=-a * ph / var)


def fmt(v, f="{:.3f}"):
    return " ".join(f.format(x) for x in v)


print(f"net {net}: MC {pre} (full / h0 / h1); slices diag | (2,2) on {ia.size} pairs | (3,1) off-diagonal; "
      f"chain dump {'chain_off%d_o.npz with inputs %s_h0' % (net, cpre) if ch is not None else 'not used'}")
summ = []
for l in range(0, L - 1):
    W = Wcol[l + 1]; H = W * W
    F = M["full"]
    t = truth(F, l); tn = nm(t)
    R = {}
    for h in ("full", "h0", "h1"):
        th = t if h == "full" else truth(M[h], l)
        R[h] = sub(th, T_pair(W, H, *ytens(M[h], l)))
    Rf = R["full"]
    rel = [a / b for a, b in zip(nm(Rf), tn)]
    sc = [float(np.corrcoef(a.ravel(), b.ravel())[0, 1]) for a, b in zip(R["h0"], R["h1"])]
    nf = [(np.linalg.norm(a - b) / 2) ** 2 / np.linalg.norm(c) ** 2 for a, b, c in zip(R["h0"], R["h1"], Rf)]
    # candidates
    vy = F["var_y"][l].astype(np.float64); Cy = F["cov_y"][l].astype(np.float64); Cy = 0.5 * (Cy + Cy.T)
    vz1 = F["var"][l + 1].astype(np.float64); Cz1 = sym0(F["cov"][l + 1])
    so2 = vz1 - H @ vy; sd2 = H @ vy
    A = (2 * so2, (so2[ia] + so2[ib]) / 3, Cz1)
    S = W @ Cy @ W.T; Sd = np.diag(S).copy()
    Bfull = (3 * Sd * Sd, Sd[ia] * Sd[ib] + 2 * S[ia, ib] ** 2, offd(3 * Sd[:, None] * S))
    Cyo = offd(Cy)
    Bsh = sub(Bfull, T_pair(W, H, 3 * vy * vy, sym0(np.outer(vy, vy) + 2 * Cyo * Cyo), offd(3 * vy[:, None] * Cyo)))
    g = gate(F["mu"][l].astype(np.float64), F["var"][l].astype(np.float64))
    Cz = sym0(F["cov"][l]); Ph, rho = g["Ph"], g["dl"]
    e2 = 2 * Ph - 2 * g["m"] * rho - 2 * Ph * Ph; f2 = 2 * g["m"] * (1 - Ph)
    G = Cz @ (W * Ph[None, :]).T                     # G_ai = sum_b C_ab Phi_b W_ib
    C2 = Cz * Cz
    Y = W * rho[None, :] * G.T                        # Y_ib = W_ib rho_b G_bi
    t1 = G.T ** 2 - (H * (Ph * Ph)[None, :]) @ C2     # (i, a): sum_(b != c) x_b x_c, x_b = W_ib Phi_b C_ab
    t2 = (Cz @ Y.T).T - W * Ph[None, :] * ((W * rho[None, :]) @ C2)
    Cd = 6 * np.einsum("ia,ia->i", H * e2[None, :], t1) + 12 * np.einsum("ia,ia->i", H * f2[None, :], t2)

    # D: the input-radius mode
    mu1 = F["mu"][l + 1].astype(np.float64); S1 = F["cov"][l + 1].astype(np.float64); S1 = 0.5 * (S1 + S1.T)
    k31 = F["k3"][l + 1].astype(np.float64); D1 = F["D21"][l + 1].astype(np.float64); np.fill_diagonal(D1, 0.0)
    my = F["mu_y"][l].astype(np.float64); k3y = F["k3_y"][l].astype(np.float64); Dy = offd(F["D21_y"][l])
    S1d = np.diag(S1).copy()
    Dfull = (GS * 3 * S1d ** 2 + GM * 4 * mu1 * k31,
             GS * (S1d[ia] * S1d[ib] + 2 * S1[ia, ib] ** 2) + GM * (2 * mu1[ia] * D1[ib, ia] + 2 * mu1[ib] * D1[ia, ib]),
             offd(GS * 3 * S1d[:, None] * S1 + GM * (3 * mu1[:, None] * D1 + np.outer(k31, mu1))))
    Dsh = sub(Dfull, T_pair(W, H, GS * 3 * vy * vy + GM * 4 * my * k3y,
                            sym0(GS * (np.outer(vy, vy) + 2 * Cyo * Cyo) + GM * (2 * my[:, None] * Dy.T + 2 * my[None, :] * Dy)),
                            offd(GS * 3 * vy[:, None] * Cyo + GM * (3 * my[:, None] * Dy + np.outer(k3y, my)))))

    # kappa3: t3 = (D3, D21 slices of z') = T3_pair(y) + R3, R3 = transport of the all-distinct y class (the legs' content)
    def T3_pair(k3v, D):
        G3 = D @ W.T                                            # G3_ai = sum_b D_ab W_ib, D_ab = kappa(y_a, y_a, y_b)
        d3 = (W * H) @ k3v + 3 * np.einsum("ia,ai->i", H, G3)
        d21 = (H * k3v[None, :]) @ W.T + H @ D @ W.T + 2 * (W * G3.T) @ W.T
        np.fill_diagonal(d21, 0.0)
        return d3, d21
    R3 = {}
    for h in ("full", "h0", "h1"):
        S_ = M[h]
        t3h = (S_["k3"][l + 1].astype(np.float64), offd(S_["D21"][l + 1]))
        R3[h] = sub(t3h, T3_pair(S_["k3_y"][l].astype(np.float64), offd(S_["D21_y"][l])))
    t3 = (F["k3"][l + 1].astype(np.float64), offd(F["D21"][l + 1]))
    rel3 = [float(np.linalg.norm(a) / np.linalg.norm(b)) for a, b in zip(R3["full"], t3)]
    sc3 = [float(np.corrcoef(a.ravel(), b.ravel())[0, 1]) for a, b in zip(R3["h0"], R3["h1"])]
    nf3 = [float((np.linalg.norm(a - b) / 2) ** 2 / np.linalg.norm(c) ** 2) for a, b, c in zip(R3["h0"], R3["h1"], R3["full"])]
    # radial-type check on kappa3: the scale tangent's kappa3 part is Cov(s, s^2) Sym3(mu, Sigma), Cov(s, s^2) ~ (1/2) Var(s^2)
    k3full = (3 * mu1 * S1d, mu1[:, None] * 0 + (S1d[:, None] * mu1[None, :] + 2 * mu1[:, None] * S1))
    k3pair = T3_pair(3 * my * vy, offd(vy[:, None] * my[None, :] + 2 * my[:, None] * Cyo))
    E3 = sub((k3full[0], offd(k3full[1])), k3pair)
    bE3 = [float(np.vdot(r, x) / np.vdot(x, x)) for r, x in zip(R3["full"], E3)]
    efE3 = [1 - float(np.linalg.norm(r - b * x) ** 2 / np.linalg.norm(r) ** 2) for r, x, b in zip(R3["full"], E3, bE3)]

    def fit(X, R_, wts=None):
        b = [float(np.vdot(r, x) / np.vdot(x, x)) for r, x in zip(R_, X)]
        ef = [1 - float(np.linalg.norm(r - bb * x) ** 2 / np.linalg.norm(r) ** 2) for r, x, bb in zip(R_, X, b)]
        w = [1 / q ** 2 for q in tn]
        bj = sum(wi * float(np.vdot(r, x)) for wi, r, x in zip(w, R_, X)) / sum(wi * float(np.vdot(x, x)) for wi, x in zip(w, X))
        efj = [1 - float(np.linalg.norm(r - bj * x) ** 2 / np.linalg.norm(r) ** 2) for r, x in zip(R_, X)]
        return b, ef, bj, efj
    bA, efA, bAj, efAj = fit(A, Rf)
    bB, efB, bBj, efBj = fit(Bsh, Rf)
    Rd = Rf[0]
    bC = float(Rd @ Cd / (Cd @ Cd)); efC = 1 - float(np.linalg.norm(Rd - bC * Cd) ** 2 / (Rd @ Rd))
    efC1 = 1 - float(np.linalg.norm(Rd - Cd) ** 2 / (Rd @ Rd))

    def fit2(X1, X2):
        Z = np.stack([X1, X2], 1); c, *_ = np.linalg.lstsq(Z, Rd, rcond=None)
        return c, 1 - float(np.linalg.norm(Rd - Z @ c) ** 2 / (Rd @ Rd))
    cAC, efAC = fit2(A[0], Cd); cBC, efBC = fit2(Bsh[0], Cd)
    bD, efD, bDj, efDj = fit(Dsh, Rf)
    efD1 = [1 - float(np.linalg.norm(r - x) ** 2 / np.linalg.norm(r) ** 2) for r, x in zip(Rf, Dsh)]
    cDC, efDC = fit2(Dsh[0], Cd)
    efDC1 = 1 - float(np.linalg.norm(Rd - Dsh[0] - Cd) ** 2 / (Rd @ Rd))
    Rres = sub(Rf, Dsh); bBr, efBr, bBrj, efBrj = fit(Bsh, Rres)
    # gain mode: Var|y|^2 from pair-supported y statistics (Cov(y_a^2, y_b^2) expanded about the means)
    k4y = F["k4_y"][l].astype(np.float64); Ky = sym0(F["K22_y"][l])
    cyy = Ky + 2 * Cyo * Cyo + 4 * np.outer(my, my) * Cyo + 2 * my[:, None] * Dy.T + 2 * my[None, :] * Dy
    var_raw = cyy.sum() + np.sum(k4y + 2 * vy * vy + 4 * my * my * vy + 4 * my * k3y)
    var_cen = (Ky + 2 * Cyo * Cyo).sum() + np.sum(k4y + 2 * vy * vy)
    g_raw = var_raw / np.sum(vy + my * my) ** 2; g_cen = var_cen / np.sum(vy) ** 2
    g_conn = (np.sum(k4y) + Ky.sum() + 4 * (my @ k3y + np.sum(Dy @ my))) / np.sum(vy + my * my) ** 2
    cG = g_conn / GS
    efG = [1 - float(np.linalg.norm(r - cG * x) ** 2 / np.linalg.norm(r) ** 2) for r, x in zip(Rf, Dsh)]
    corrD = [float(np.corrcoef(r.ravel(), x.ravel())[0, 1]) for r, x in zip(Rf, Dsh)]
    line = (f"layer {l:2d}->{l + 1:2d}: |R|/|t| {fmt(rel)}  split corr {fmt(sc, '{:+.2f}')}  noise share {fmt(nf)}  "
            f"|Cgauss|/|t_d| {np.linalg.norm(Cd) / tn[0]:.3f}\n"
            f"    A lam-core : lam* {fmt(bA, '{:.5f}')} joint {bAj:.5f}; explained {fmt(efA)} (joint {fmt(efAj)})\n"
            f"    B scale    : g {fmt(bB, '{:.4f}')} joint {bBj:.4f}; explained {fmt(efB)} (joint {fmt(efBj)}); "
            f"lam*/(3 g <sd2>) {bAj / (3 * bBj * sd2.mean()) if bBj else float('nan'):.3f}\n"
            f"    C gauss2 (diag): coef {bC:+.3f} explained {efC:.3f} (coef fixed 1: {efC1:+.3f}); "
            f"A+C {efAC:.3f} (coefs {cAC[0]:.5f} {cAC[1]:+.3f}); B+C {efBC:.3f} (coefs {cBC[0]:.4f} {cBC[1]:+.3f})\n"
            f"    D radius   : |D|/|t| {fmt([a / b for a, b in zip(nm(Dsh), tn)])} corr {fmt(corrD, '{:+.3f}')}; coef 1: explained {fmt(efD1, '{:+.3f}')}; "
            f"fitted {fmt(bD)} joint {bDj:.3f} (explained {fmt(efDj, '{:+.3f}')}); D+C diag {efDC:.3f} (coefs {cDC[0]:+.3f} {cDC[1]:+.3f}; "
            f"both fixed 1: {efDC1:+.3f}); after D, scale-mode g {bBrj:.5f} (2/n = {GS:.5f}) explains {fmt(efBrj, '{:+.3f}')}\n"
            f"    gain mode: g_raw {g_raw:.5f} (-2/n: {g_raw - GS:.5f}), g_cen {g_cen:.5f} (-2/n: {g_cen - GS:.5f}) | fitted g_B {bBj:.5f}, "
            f"D coef x Var(r^2) {bDj * GS:.5f}; g_conn {g_conn:.5f} -> D coef {cG:.3f}, explained {fmt(efG, '{:+.3f}')}")
    line += (f"\n    kappa3: |R3|/|t3| (D3 | D21) {fmt(rel3)} split corr {fmt(sc3, '{:+.2f}')} noise share {fmt(nf3)}; "
             f"scale-mode shape fit coef (x2 = Var(s^2)) {fmt(bE3, '{:.4f}')} explained {fmt(efE3, '{:+.3f}')}")
    rec = dict(l=l, rel=rel, sc=sc, nf=nf, rel3=rel3, sc3=sc3, efE3=efE3, bAj=bAj, efAj=efAj, bBj=bBj, efBj=efBj, efC=efC, bC=bC, efAC=efAC, efBC=efBC,
               efD1=efD1, bDj=bDj, efDC1=efDC1, bBrj=bBrj, g_raw=g_raw, g_cen=g_cen,
               g_conn=g_conn, efG=efG, efDj=efDj)
    need = [f"K22_{l}", f"K4v_{l}", f"K2v_{l}", f"var_{l}", f"C_off_{l}", f"D21_{l}", f"wk4m_{l}", f"wk431_{l}",
            f"g4rowown_{l + 1}", f"wk4mown_{l + 1}", f"wk431own_{l + 1}", f"C_off_{l + 1}", f"var_{l + 1}"]
    if ch is not None and all(k in ch.files for k in need):
        d_c = ch[f"K4v_{l}"].astype(np.float64); K_c = sym0(ch[f"K22_{l}"])
        dY, KY, BY = ytens(F, l)
        gc = gate(I0["mu"][l].astype(np.float64), ch[f"var_{l}"].astype(np.float64))
        C0 = offd(ch[f"C_off_{l}"]); Gm = offd(ch[f"D21_{l}"]); Kz = offd(ch[f"wk4m_{l}"]); Bz = offd(ch[f"wk431_{l}"].astype(np.float64).T)
        Pc, dl, dl1, m, v = gc["Ph"], gc["dl"], gc["dl1"], gc["m"], gc["v"]
        e31 = 3 * (v - m * m) * (1 - Pc); e32 = 6 * m * (1 - Pc) + 3 * (m * m - v) * dl
        e33 = 6 * Pc - 6 * m * dl + 3 * m * m * dl1 - 3 * v * dl1
        By = offd(e31[:, None] * Pc[None, :] * C0 + 0.5 * Gm * e32[:, None] * Pc[None, :] + 0.5 * Gm.T * e31[:, None] * dl[None, :]
                  + 0.25 * Kz * e32[:, None] * dl[None, :] + Bz / 6 * e33[:, None] * Pc[None, :] + Bz.T / 6 * e31[:, None] * dl1[None, :])
        yrel = lambda p, q: (float(np.linalg.norm(p - q) / np.linalg.norm(q)), float(np.corrcoef(p.ravel(), q.ravel())[0, 1]))
        yd, yK, yB = yrel(d_c, dY), yrel(K_c[~np.eye(n, dtype=bool)], KY[~np.eye(n, dtype=bool)]), yrel(By[~np.eye(n, dtype=bool)], BY[~np.eye(n, dtype=bool)])
        ed = T_pair(W, H, d=d_c - dY); eK = T_pair(W, H, K=K_c - KY); eB = T_pair(W, H, B=By - BY)
        Cz_c = ch[f"C_off_{l + 1}"].astype(np.float64); vz_c = ch[f"var_{l + 1}"].astype(np.float64)
        own431 = ch[f"wk431own_{l + 1}"].astype(np.float64); offm = ~np.eye(n, dtype=bool)
        lam = float(np.sum(own431[offm] * Cz_c[offm]) / np.sum(Cz_c[offm] ** 2))
        so2c = vz_c - H @ ch[f"K2v_{l}"].astype(np.float64)
        lcore = (2 * lam * so2c, lam * (so2c[ia] + so2c[ib]) / 3, lam * offd(Cz_c))
        own = (ch[f"g4rowown_{l + 1}"].astype(np.float64), ch[f"wk4mown_{l + 1}"].astype(np.float64)[ia, ib], offd(own431.T))
        c6x = T_pair(W, H, d_c, K_c, By)
        clos = sub(c6x, T_pair(W, H, dY, KY, BY)); cls = sub(lcore, Rf); rep = sub(own, add(c6x, lcore))
        ech = sub(own, t)
        cls_g = sub(tuple(bBj * x for x in Bsh), Rf); cls_D = sub(Dsh, Rf)
        # compensation: cosine between the closure error and the lam core's misfit of the omitted classes (note XXI s7)
        cosCC = [float(np.vdot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))) for a, b in zip(clos, cls)]
        ech_G = add(clos, sub(tuple(cG * x for x in Dsh), Rf), rep)        # lam core replaced by g_conn D
        ech_Gf = add(clos, sub(tuple(bDj * x for x in Dsh), Rf), rep)      # ... by the best-fit coefficient of D
        ech_0 = add(clos, sub((0 * Rf[0], 0 * Rf[1], 0 * Rf[2]), Rf), rep)  # ... by nothing
        r_ = lambda X: fmt([a / b for a, b in zip(nm(X), tn)])
        line += (f"\n    y slices (chain from true z vs MC): K4v {yd[0]:.3f}/{yd[1]:+.3f}  K22 {yK[0]:.3f}/{yK[1]:+.3f}  "
                 f"B^y closed form {yB[0]:.3f}/{yB[1]:+.3f}; transported: K4v {r_(ed)} | K22 {r_(eK)} | B^y {r_(eB)}\n"
                 f"    chain error {r_(ech)} = closure {r_(clos)} + classes-lam {r_(cls)} + repr {r_(rep)}; "
                 f"chain lam {lam:.5f}; classes with best scale mode {r_(cls_g)}; classes with radius mode D (no lam) {r_(cls_D)}\n"
                 f"    compensation cos(closure, lam core - R) {fmt(cosCC, '{:+.2f}')}; chain error with the lam core replaced by: "
                 f"g_conn D {r_(ech_G)} | best-fit D {r_(ech_Gf)} | nothing {r_(ech_0)}")
        rec.update(lam=lam, ech=[a / b for a, b in zip(nm(ech), tn)], clos=[a / b for a, b in zip(nm(clos), tn)],
                   cls=[a / b for a, b in zip(nm(cls), tn)], cosCC=cosCC, echG=[a / b for a, b in zip(nm(ech_G), tn)],
                   echGf=[a / b for a, b in zip(nm(ech_Gf), tn)], ech0=[a / b for a, b in zip(nm(ech_0), tn)])
    print(line, flush=True)
    summ.append(rec)
print("summary (mean over transitions 2->3 .. 14->15):")
S_ = [r for r in summ if r["l"] >= 2]
for key in ("rel3", "sc3", "efE3", "rel", "sc", "nf", "efAj", "efBj", "efD1", "efDj", "efG", "ech", "clos", "cls", "cosCC", "echG", "echGf", "ech0"):
    if all(key in r for r in S_):
        print(f"  {key:5s} " + fmt(np.mean([r[key] for r in S_], 0)))
for key in ("bAj", "bBj", "g_raw", "g_cen", "g_conn", "bDj", "bBrj", "efC", "bC", "efAC", "efBC", "efDC1", "lam"):
    if all(key in r for r in S_):
        print(f"  {key:5s} " + " ".join(f"{r[key]:.4g}" for r in S_))
