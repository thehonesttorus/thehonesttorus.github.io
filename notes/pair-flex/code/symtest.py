# Test S of note XXXVIII (notes/pair-flex): is the transport R of the omitted post-activation classes ((2,1,1) and
# (1,1,1,1)) equal to the transport of their own collective projections?
#   python symtest.py NET "SYMCHUNK_GLOB" [KLIST] [MCPREFIX]          e.g. python symtest.py 1 "sym1_*.npz" 8,16,32,64 mc4_off1
# SYMCHUNKs from mcsym.py (half 0 = even chunk ids, half 1 = odd). R per transition l -> l+1 from MCPREFIX_{full,h0,h1}.npz
# (mcstats.py ... post): R = t - T_pair(y), exactly W# of the y classes the pair state omits (note XXXVI section 3).
# From the raw sums, with every index coincidence removed by inclusion-exclusion and every mean corrected exactly:
#   N_a^off = U^T Psi_a U, Psi_a = unit a's (2,1,1) block kappa(e_a, e_a, e_c, e_d) (c != d, both != a), K x K;
#   T4^off  = the all-distinct class contracted with U on all four indices, K^4.
# The model tensor Yhat = T4^off[U^4] + sum_(6 slot pairs) delta_pair (U M_a U^T), M_a = N_a^off - T4^off[U_a, U_a], equals
# the collective projections on the omitted classes; R_pred = W#[Yhat] - T_pair(Yhat's own pair slices), parameter-free.
import sys, os, glob, re, math, numpy as np


def load_chunks(files):
    tot = None
    for f in files:
        z = np.load(f)
        if tot is None:
            tot = {k: (np.array(z[k], np.float64) if z[k].dtype.kind == "f" and not k.startswith(("U_", "my_")) else np.array(z[k]))
                   for k in z.files}
        else:
            for k in z.files:
                if k.startswith(("U_", "my_")) or k in ("K", "layers"):
                    continue
                tot[k] = tot[k] + z[k]
    return tot


def pair_index(K):
    iu = np.triu_indices(K); PI = np.zeros((K, K), int)
    PI[iu] = np.arange(len(iu[0])); PI[(iu[1], iu[0])] = np.arange(len(iu[0]))
    return iu, PI


def symbols(acc, l, K):
    """N^off (n x K x K), T4^off (K^4), U (n x K), plus pieces for diagnostics, from mcsym sums at layer l, first K modes."""
    Kmax = int(acc["K"]); c = float(acc["c"])
    iuF = np.triu_indices(Kmax); selP = np.where((iuF[0] < K) & (iuF[1] < K))[0]
    g = lambda k: acc[f"{k}_{l}"]
    U = np.asarray(acc[f"U_{l}"], np.float64)[:, :K]
    s1, s2, s3, s4 = (g(k) / c for k in ("s1", "s2", "s3", "s4"))
    Craw, Draw = g("C") / c, g("D") / c                       # Draw[c, a] = E[e_c^2 e_a]
    et, e2t, e3t = g("et")[:, :K] / c, g("e2t")[:, :K] / c, g("e3t")[:, :K] / c
    etp, e2tp, e2S = g("etp")[:, selP] / c, g("e2tp")[:, selP] / c, g("e2S")[:, selP] / c
    t1, tp = g("t1")[:K] / c, g("tp")[selP] / c
    tpt = g("tpt")[selP][:, :K] / c; tptp = g("tptp")[np.ix_(selP, selP)] / c
    S1, SS = g("S1")[selP] / c, g("SS")[np.ix_(selP, selP)] / c
    iu, PI = pair_index(K); k_, m_ = iu
    UU = U[:, k_] * U[:, m_]                                  # n x q
    eb, tb = s1, t1
    Ch = Craw - np.outer(eb, eb); var = np.diag(Ch).copy()
    mu3 = s3 - 3 * eb * s2 + 2 * eb ** 3
    mu4 = s4 - 4 * eb * s3 + 6 * eb ** 2 * s2 - 3 * eb ** 4
    k4 = mu4 - 3 * var ** 2
    Xt = Ch @ U                                               # E[x_a tau_k]
    Lam = U.T @ Ch @ U
    # N^full_a = kappa(x_a, x_a, tau_k, tau_m), pair form (n x q)
    tk, tm = tb[k_][None, :], tb[m_][None, :]
    e_ = eb[:, None]
    mu22 = (e2tp - 2 * e_ * etp - tk * e2t[:, m_] - tm * e2t[:, k_] + e_ ** 2 * tp[None, :]
            + 2 * e_ * tk * et[:, m_] + 2 * e_ * tm * et[:, k_] + s2[:, None] * tk * tm - 3 * e_ ** 2 * tk * tm)
    Nfull = mu22 - var[:, None] * Lam[k_, m_][None, :] - 2 * Xt[:, k_] * Xt[:, m_]
    # sum_c Khat_ac U_ck U_cm, Khat_ac = kappa(x_a, x_a, x_c, x_c) (c = a included)
    eUU = eb[:, None] * UU
    s0 = UU.T @ (eb * eb)
    ExS = (e2S - 2 * e_ * (Draw.T @ UU) + e_ ** 2 * S1[None, :] - 2 * (Draw @ eUU) + 4 * e_ * (Craw @ eUU)
           - 2 * e_ ** 2 * s0[None, :] + var[:, None] * s0[None, :])
    Khs = ExS - var[:, None] * (UU.T @ var)[None, :] - 2 * (Ch * Ch) @ UU
    # g~_a = sum_(d != a) kappa(x_a, x_a, x_a, x_d) U_d
    mu31 = (e3t - 3 * e_ * e2t + 3 * e_ ** 2 * et - 3 * e_ ** 3 * tb[None, :] - tb[None, :] * s3[:, None]
            + 3 * e_ * tb[None, :] * s2[:, None])
    k31t = mu31 - 3 * var[:, None] * Xt
    gt = k31t - k4[:, None] * U
    Noff_p = Nfull - Khs - (U[:, k_] * gt[:, m_] + gt[:, k_] * U[:, m_])
    Noff = Noff_p[:, PI]                                      # n x K x K
    # T4^off = kappa4(tau) - [one pair] - [two pairs] - [triple] - [all four]
    M2 = tp[PI]; M3 = tpt[PI]; M4 = tptp[PI[:, :, None, None], PI[None, None, :, :]]
    ein = np.einsum
    mu4t = (M4 - ein("k,mpq->kmpq", tb, M3) - ein("m,kpq->kmpq", tb, M3) - ein("p,kmq->kmpq", tb, M3) - ein("q,kmp->kmpq", tb, M3)
            + ein("k,m,pq->kmpq", tb, tb, M2) + ein("k,p,mq->kmpq", tb, tb, M2) + ein("k,q,mp->kmpq", tb, tb, M2)
            + ein("m,p,kq->kmpq", tb, tb, M2) + ein("m,q,kp->kmpq", tb, tb, M2) + ein("p,q,km->kmpq", tb, tb, M2)
            - 3 * ein("k,m,p,q->kmpq", tb, tb, tb, tb))
    k4t = mu4t - (ein("km,pq->kmpq", Lam, Lam) + ein("kp,mq->kmpq", Lam, Lam) + ein("kq,mp->kmpq", Lam, Lam))
    N4 = (UU.T @ Noff_p)[PI[:, :, None, None], PI[None, None, :, :]]          # sum_a U_ak U_am Noff_a[p, q]
    one_pair = (N4 + ein("kpmq->kmpq", N4) + ein("kqmp->kmpq", N4) + ein("mpkq->kmpq", N4) + ein("mqkp->kmpq", N4)
                + ein("pqkm->kmpq", N4))
    ESx = SS - 2 * (UU.T @ Draw @ eUU) - 2 * (UU.T @ Draw @ eUU).T + 4 * (eUU.T @ Craw @ eUU) + np.outer(S1, s0) + np.outer(s0, S1) - 3 * np.outer(s0, s0)
    Sxm = S1 - s0
    UKU = ESx - np.outer(Sxm, Sxm) - 2 * (UU.T @ (Ch * Ch) @ UU)
    F4p = UU.T @ (k4[:, None] * UU)
    TP = (UKU - F4p)[PI[:, :, None, None], PI[None, None, :, :]]
    two_pair = TP + ein("kpmq->kmpq", TP) + ein("kqmp->kmpq", TP)
    G3 = (UU.T @ (U[:, :, None] * gt[:, None, :]).reshape(len(U), -1)).reshape(-1, K, K)[PI]   # [k,m,p,q] = sum_a U_ak U_am U_ap gt_aq
    triple = G3 + ein("kmqp->kmpq", G3) + ein("kpqm->kmpq", G3) + ein("mpqk->kmpq", G3)
    F4 = F4p[PI[:, :, None, None], PI[None, None, :, :]]
    T4off = k4t - one_pair - two_pair - triple - F4
    return dict(U=U, Noff=Noff, T4off=T4off, Nfull=Nfull[:, PI], k4t=k4t, var=var, eb=eb, Lam=Lam)


def T_pair(W, H, ia, ib, d=None, K=None, B=None):
    """exact slices of W#[D(d) + S22(K) + S31(B)] (yclasses.py): diag, (2,2) on pairs (ia, ib), (3,1)_ij (zero diagonal)."""
    n = W.shape[0]
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


def wsharp_sym(W, H, ia, ib, U, M, cross=True):
    """W# of the (2,1,1) entries of sum_(6 slot pairs) delta_pair (U M_a U^T): slices diag, (2,2) on pairs, (3,1)."""
    n, K = U.shape
    A = W @ U; Mf = M.reshape(n, K * K)
    AA = (A[:, :, None] * A[:, None, :]).reshape(n, K * K)
    Q = AA @ Mf.T                                              # Q_ia = A_i^T M_a A_i
    MM = (H @ Mf).reshape(n, K, K)                             # sum_a W_ia^2 M_a
    V = np.einsum("ikm,im->ik", MM, A)
    diag = 6 * np.sum(H * Q, axis=1)
    k31 = 3 * (V @ A.T + (W * Q) @ W.T)
    d22 = np.einsum("pk,pkm,pm->p", A[ib], MM[ia], A[ib]) + np.einsum("pk,pkm,pm->p", A[ia], MM[ib], A[ia])
    x22 = np.zeros(ia.size)
    if cross:
        for s in range(0, ia.size, 2000):
            sl = slice(s, s + 2000)
            Ms = ((W[ia[sl]] * W[ib[sl]]) @ Mf).reshape(-1, K, K)
            x22[sl] = 4 * np.einsum("pk,pkm,pm->p", A[ia[sl]], Ms, A[ib[sl]])
    k22 = d22 + x22
    # own pair slices of the symbol tensor, removed: diag 6 U_a^T M_a U_a, (2,2) U_b^T M_a U_b + U_a^T M_b U_a, (3,1) 3 U_a^T M_a U_b
    UUf = (U[:, :, None] * U[:, None, :]).reshape(n, K * K)
    qab = UUf @ Mf.T                                           # qab[b, a] = U_b^T M_a U_b
    od = 6 * np.einsum("ak,akm,am->a", U, M, U)
    oK = qab + qab.T; np.fill_diagonal(oK, 0.0)
    oB = 3 * np.einsum("akm,am->ak", M, U) @ U.T; np.fill_diagonal(oB, 0.0)
    pd, p22, p31 = T_pair(W, H, ia, ib, od, oK, oB)
    np.fill_diagonal(k31, 0.0)
    return (diag - pd, k22 - p22, k31 - p31), (np.linalg.norm(x22), np.linalg.norm(k22))


def wsharp_t4(W, H, ia, ib, U, T):
    """W# of the (2,1,1) and (1,1,1,1) entries of T[U^4]."""
    n, K = U.shape
    def slices(A):
        X = ((A[:, :, None] * A[:, None, :]).reshape(n, K * K) @ T.reshape(K * K, K * K)).reshape(n, K, K)   # T[A_i, A_i, ., .]
        dg = np.einsum("ikm,ik,im->i", X, A, A)
        Y = np.einsum("ikm,im->ik", X, A); b31 = Y @ A.T; np.fill_diagonal(b31, 0.0)
        return X, dg, b31
    A = W @ U
    X, dg, b31 = slices(A)
    k22 = np.einsum("pkm,pk,pm->p", X[ia], A[ib], A[ib])
    Xu, ou, Bu = slices(U)
    Ku = np.einsum("akm,bk,bm->ab", Xu, U, U); np.fill_diagonal(Ku, 0.0)
    pd, p22, p31 = T_pair(W, H, ia, ib, ou, Ku, Bu)
    return dg - pd, k22 - p22, b31 - p31


def r_pred(W, H, ia, ib, S, cross=True):
    U, Noff, T4 = S["U"], S["Noff"], S["T4off"]
    T4UU = np.einsum("kmpq,ak,am->apq", T4, U, U)
    r211, cr = wsharp_sym(W, H, ia, ib, U, Noff, cross)
    r_t4 = wsharp_t4(W, H, ia, ib, U, T4)
    r_t4_211, _ = wsharp_sym(W, H, ia, ib, U, T4UU, cross)
    r1111 = tuple(a - b for a, b in zip(r_t4, r_t4_211))
    return tuple(a + b for a, b in zip(r211, r1111)), r211, r1111, cr


if __name__ == "__main__":
    import math as _m
    net, pat = int(sys.argv[1]), sys.argv[2]
    KL = [int(x) for x in (sys.argv[3] if len(sys.argv) > 3 else "8,16,32,64").split(",")]
    pre = sys.argv[4] if len(sys.argv) > 4 else f"mc4_off{net}"
    files = sorted(glob.glob(pat)); idx = {int(re.search(r"_(\d+)\.npz$", f).group(1)): f for f in files}
    sym = {"h0": load_chunks([f for i, f in idx.items() if i % 2 == 0]), "h1": load_chunks([f for i, f in idx.items() if i % 2 == 1])}
    sym["full"] = {k: (sym["h0"][k] + sym["h1"][k] if not (k.startswith(("U_", "my_")) or k in ("K", "layers")) else sym["h0"][k])
                   for k in sym["h0"]}
    Wcol = np.load(f"../official/W_off{net}.npy").astype(np.float64); L, n, _ = Wcol.shape
    Mc = {h: np.load(f"{pre}_{h}.npz") for h in ("full", "h0", "h1")}
    rng = np.random.default_rng(1)
    ia = rng.integers(0, n, 12000); ib = rng.integers(0, n, 12000); k_ = ia != ib; ia, ib = ia[k_], ib[k_]
    offd = lambda X: (lambda Y: (np.fill_diagonal(Y, 0.0), Y)[1])(np.array(X, np.float64))
    sym0 = lambda X: (lambda Y: (np.fill_diagonal(Y, 0.0), Y)[1])(0.5 * (np.array(X, np.float64) + np.array(X, np.float64).T))
    en = lambda X: float(np.sum(np.asarray(X) ** 2))
    print(f"net {net}: symbols from {len(files)} chunks ({int(sym['full']['c'])} inputs), R from {pre}; slices diag | (2,2) on "
          f"{ia.size} pairs | (3,1) off-diagonal")
    GS = 2.0 / n; Er = _m.sqrt(2.0 / n) * _m.exp(_m.lgamma((n + 1) / 2) - _m.lgamma(n / 2)); GM = (1 + 2.0 / n) - Er * Er * (n + 1.0) / n
    for l in [int(x) for x in sym["full"]["layers"]]:
        W = Wcol[l + 1]; H = W * W
        R = {}
        for h in Mc:
            F = Mc[h]
            t = (F["k4"][l + 1].astype(np.float64), F["K22"][l + 1].astype(np.float64)[ia, ib], offd(F["K31"][l + 1]))
            y = (F["k4_y"][l].astype(np.float64), sym0(F["K22_y"][l]), offd(F["K31_y"][l]))
            R[h] = tuple(a - b for a, b in zip(t, T_pair(W, H, ia, ib, *y)))
        nR = [en(a - b) / 4 for a, b in zip(R["h0"], R["h1"])]
        # candidate D (yclasses.py), fitted per slice
        F = Mc["full"]
        vy = F["var_y"][l].astype(np.float64); Cy = F["cov_y"][l].astype(np.float64); Cyo = offd(0.5 * (Cy + Cy.T))
        mu1 = F["mu"][l + 1].astype(np.float64); S1 = F["cov"][l + 1].astype(np.float64); S1 = 0.5 * (S1 + S1.T)
        k31 = F["k3"][l + 1].astype(np.float64); D1 = F["D21"][l + 1].astype(np.float64); np.fill_diagonal(D1, 0.0)
        my = F["mu_y"][l].astype(np.float64); k3y = F["k3_y"][l].astype(np.float64); Dy = offd(F["D21_y"][l]); S1d = np.diag(S1).copy()
        Dfull = (GS * 3 * S1d ** 2 + GM * 4 * mu1 * k31,
                 GS * (S1d[ia] * S1d[ib] + 2 * S1[ia, ib] ** 2) + GM * (2 * mu1[ia] * D1[ib, ia] + 2 * mu1[ib] * D1[ia, ib]),
                 offd(GS * 3 * S1d[:, None] * S1 + GM * (3 * mu1[:, None] * D1 + np.outer(k31, mu1))))
        Dsh = tuple(a - b for a, b in zip(Dfull, T_pair(W, H, ia, ib, GS * 3 * vy * vy + GM * 4 * my * k3y,
                    sym0(GS * (np.outer(vy, vy) + 2 * Cyo * Cyo) + GM * (2 * my[:, None] * Dy.T + 2 * my[None, :] * Dy)),
                    offd(GS * 3 * vy[:, None] * Cyo + GM * (3 * my[:, None] * Dy + np.outer(k3y, my))))))
        Rf = R["full"]
        def ex(X, nX=None):           # explained share at coefficient 1, raw and noise-corrected (nX: noise energy of X)
            raw = [1 - en(r - x) / en(r) for r, x in zip(Rf, X)]
            if nX is None:
                return raw, None
            cor = [1 - (en(r - x) - a - b) / (en(r) - a) for r, x, a, b in zip(Rf, X, nR, nX)]
            return raw, cor
        def fitc(X):
            return [float(np.vdot(r, x) / np.vdot(x, x)) for r, x in zip(Rf, X)]
        bD = fitc(Dsh); exD = [1 - en(r - b * x) / en(r) for r, x, b in zip(Rf, Dsh, bD)]
        save = dict(R0=Rf[0], R1=Rf[1], R2=Rf[2], D0=Dsh[0], D1=Dsh[1], D2=Dsh[2], ia=ia, ib=ib)
        f3 = lambda v, f="{:+.3f}": " ".join(f.format(x) for x in v)
        print(f"\nlayer {l:2d} -> {l + 1:2d}: |R| {f3([np.sqrt(en(r)) for r in Rf], '{:.3e}')}; R noise share {f3([a / en(r) for a, r in zip(nR, Rf)], '{:.3f}')}; "
              f"candidate D fitted: coef {f3(bD, '{:.2f}')} explains {f3(exD)}", flush=True)
        for K in KL:
            Sf = symbols(sym["full"], l, K); Sh = [symbols(sym[h], l, K) for h in ("h0", "h1")]
            P, p211, p1111, cr = r_pred(W, H, ia, ib, Sf)
            save.update({f"P{K}_{j}": P[j] for j in range(3)})
            Ph = [r_pred(W, H, ia, ib, s)[0] for s in Sh]
            nP = [en(a - b) / 4 for a, b in zip(Ph[0], Ph[1])]
            raw, cor = ex(P, nP)
            b = fitc(P); corr = [float(np.corrcoef(r.ravel(), x.ravel())[0, 1]) for r, x in zip(Rf, P)]
            e211 = [1 - en(r - x) / en(r) for r, x in zip(Rf, p211)]
            ratio = [np.sqrt(en(a) / max(en(c_), 1e-300)) for a, c_ in zip(p1111, p211)]
            # jointly with D (two coefficients per slice)
            jD = []
            for r, x, dd in zip(Rf, P, Dsh):
                Z = np.stack([x.ravel(), dd.ravel()], 1); cc, *_ = np.linalg.lstsq(Z, r.ravel(), rcond=None)
                jD.append(1 - en(r.ravel() - Z @ cc) / en(r))
            print(f"  K={K:2d}: R_pred explains {f3(raw)} (noise-corrected {f3(cor)}); coef {f3(b, '{:.2f}')}; corr {f3(corr)}; "
                  f"R_pred noise share {f3([a / en(r) for a, r in zip(nP, Rf)], '{:.3f}')}\n"
                  f"         (2,1,1) part alone {f3(e211)}; |(1,1,1,1) part| / |(2,1,1) part| {f3(ratio, '{:.2f}')}; "
                  f"with D jointly {f3(jD)} (D alone {f3(exD)}); (2,2) incoherent cross term |x|/|R_pred22| {cr[0] / max(cr[1], 1e-300):.3f}",
                  flush=True)
        if os.environ.get("SYMSAVE"):
            np.savez(f"symR_off{net}_{l}.npz", **save)
