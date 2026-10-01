#!/usr/bin/env python3
"""Transfer-operator / Lyapunov-spectrum measurement of the transported third cumulant of a random ReLU MLP (numpy only).

The question (bridges stream, 2026-10-01): expander / transfer-operator intuition predicts that products of the layer maps
that carry "old content" of the third cumulant mix fast and concentrate it on a few top directions (the top of the Lyapunov
spectrum of products of random masked matrices), which would make cheap carriers possible. This script measures it.

Conventions (x@W, as in the atlases): z_l = a_{l-1} @ W[l] (a_{-1} = x ~ N(0, I)), a_l = relu(z_l); column form A_l = W[l].T,
gates Phi_l = gate_p[l] = P(z_l > 0), D_l = diag(Phi_l). The linear ("old content", births ignored) propagator from a_s to z_t
(s < t; t - s weight matrices, t - s - 1 gates) is

    J_{s->t} = A_t D_{t-1} A_{t-1} ... D_{s+1} A_{s+1}.

The gated version G_{s->t} = D_t J_{s->t} (a_s -> linearised a_t) is what streams/old-content/propproj.py calls M_{s->t}; its PR is
printed too for cross-checking. For the third cumulant the layer map on symmetric 3-tensors is X -> A_{l+1}^{(x)3}(Phi_l^{(x)3} . X),
so the transported tensor is X_{s,t} = J_{s->t}^{(x)3} K_s; for covariance-like content it is S -> B S B^T with B_l = A_{l+1} D_l
(z_l -> z_{l+1}), i.e. x@W form S -> W^T diag(Phi) S diag(Phi) W.

Parts (mode `atlas A.npz [B.npz]`; B = the same MLP with an independent sample seed, used for noise correction):
 (a) singular values of J_{s->t} for all s < t: participation ratio PR = (sum s^2)^2 / sum s^4, r90 / r99 / r999 (ranks holding
     90 / 99 / 99.9 % of sum s^2), the free-probability prediction PR_free = n / r_pred with
     r_pred - 1 = sum_l (r(A_l^T A_l) - 1) + sum_gates (r(D_l^2) - 1),  r(x) = phi(x^2) / phi(x)^2,
     finite-time Lyapunov exponents of J_{0->L-1}, and alignment of top-k singular subspaces across sources (U side) and targets
     (V side).
 (b) K_s = all-distinct part of kappa3(a_s) (central moments from post_M3, as in experiments/oracle_k3.py), transported by
     J_{s->t}^{(x)3}; its D21 contribution d_{s,t} = X_{s,t}[a, a, b] at layer t, compared with the atlas's D21(z_t). The tensor is
     projected in all three indices on the top-k left singular vectors U_k of J_{s->t} (equivalently, the source on the top-k
     right singular vectors V_k: U_k U_k^T J = J V_k V_k^T) and the relative error of the D21 contribution is reported:
        eps_own(k) = ||D21(P_k X) - d|| / ||d||,   eps_tot(k) = ||D21(P_k X) - d|| / ||D21(z_t)||,
     both noise-corrected with the independent atlas B: eps^2 = <e_A, e_B> / <d_A, d_B> (the projection is linear, the noise of the
     two atlases is independent, so the cross product removes the noise energy from numerator and denominator).
     Baselines: a random k-subspace, the data-adapted HOSVD subspace of X_{s,t} itself, and the carrier that a chain could
     actually run: the source core fixed at birth in the basis V_k(J_{s->L-1}) and transported as n x k legs.
 (c) the symmetric-matrix transfer map S -> B_l S B_l^T: its singular values on Sym_n are {s_i s_j}_{i<=j} and its eigenvalues
     {l_i l_j}_{i<=j} (s, l of B_l; the Sym^2 functor; checked numerically in `selftest`); per layer and for products; plus the
     capture of the actual transported covariance J C(a_s) J^T by U_k.
 (d) the Gaussian-ensemble prediction of eps_own(k) from the propagator alone. If the source K is replaced by a symmetric Gaussian
     tensor (or by a Haar rotation of any traceless source), E||D21 error||^2 / E||D21||^2 = 1 - f(k) / f(n) with
        f(k) = (1/3) (sum_a g_a^2)(sum_a g_a) + (2/3) sum_a g_a h_a,  g_a = sum_{p<=k} s_p^2 U_ap^2,  h_a = sum_{p<=k} s_p^4 U_ap^2.
     This needs only the SVD of J, so it runs at width 1024 (mode `scaling`), and its agreement with the measured eps_own at width
     128 is part of the atlas output.

Modes:
    python transfer_spectrum.py atlas A.npz [B.npz] [--save out.npz]
    python transfer_spectrum.py scaling [--widths 128,256,512,1024] [--depth 16] [--mc 20000] [--save out.npz]
    python transfer_spectrum.py selftest
"""
import sys
import os
import argparse
import time

sys.dont_write_bytecode = True
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "experiments"))
from oracle_k3 import central3, all_distinct, slices  # noqa: E402

KGRID = [1, 2, 3, 4, 5, 6, 8, 10, 12, 14, 16, 20, 24, 28, 32, 40, 48, 56, 64, 80, 96, 112, 128]
KCOARSE = [4, 8, 12, 16, 20, 24, 32, 40, 48, 56, 64, 80, 96, 112, 128]


# ----------------------------------------------------------------------------------------------------------- tensor helpers
def rot3(X, U):
    """Y_pqr = sum_ijk U_ip U_jq U_kr X_ijk for X (n, n, n) and U (n, m)."""
    n = X.shape[0]
    m = U.shape[1]
    T = (X.reshape(n * n, n) @ U).reshape(n, n, m)                       # [i, j, r]
    T = np.ascontiguousarray(T.transpose(0, 2, 1)).reshape(n * m, n) @ U  # [(i, r), q]
    T = (U.T @ T.reshape(n, m * m)).reshape(m, m, m)                     # [p, r, q]
    return np.ascontiguousarray(T.transpose(0, 2, 1))                    # [p, q, r]


def d21_from_core(Y, U):
    """D21 of U^{(x)3} Y for a core Y (k, k, k) and any U (n, k): d_ab = sum_pqr U_ap U_aq U_br Y_pqr."""
    n, k = U.shape
    T1 = np.ascontiguousarray(Y).reshape(k * k, k) @ U.T                  # [(p, q), b]
    UU = (U[:, :, None] * U[:, None, :]).reshape(n, k * k)
    return UU @ T1


def d21(X):
    return slices(X)[1]


def phi3(X, g):
    return X * g[:, None, None] * g[None, :, None] * g[None, None, :]


def spec_stats(sv):
    e = np.asarray(sv, dtype=np.float64) ** 2
    tot = e.sum()
    c = np.cumsum(e) / tot
    pr = tot ** 2 / (e ** 2).sum()
    rk = [int(min(np.searchsorted(c, q) + 1, len(e))) for q in (0.9, 0.99, 0.999)]
    return pr, rk[0], rk[1], rk[2]


def rmom(M):
    """r(M^T M) = phi((M^T M)^2) / phi(M^T M)^2 with phi = Tr / n."""
    n = M.shape[1]
    G = M.T @ M
    return (np.sum(G * G) / n) / (np.trace(G) / n) ** 2


def eps_ensemble(U, sv):
    """Gaussian-ensemble relative D21 error of projecting on the top-k left singular vectors, for every k = 1..n (see (d))."""
    U2 = U ** 2
    g = np.cumsum(U2 * (sv ** 2)[None, :], axis=1)        # g[a, k-1] = sum_{p<=k} s_p^2 U_ap^2
    h = np.cumsum(U2 * (sv ** 4)[None, :], axis=1)
    f = (np.sum(g ** 2, 0) * np.sum(g, 0)) / 3.0 + 2.0 * np.sum(g * h, 0) / 3.0
    return np.sqrt(np.clip(1.0 - f / f[-1], 0.0, None))


def gate_si(GG, p, lo=0.02):
    """unpinned spectral-independence constant of the gate (face) law on the uncertain gates lo < p < 1 - lo:
    eta_0 = lambda_max(Psi) with Psi(i, j) = P(g_j | g_i) - P(g_j | not g_i) = Cov_ij / Var_i (i != j), i.e. lambda_max(Cor) - 1
    (Anari-Liu-Oveis Gharan Def. 1.1; Chen-Eldan Fact 23); returns (number of uncertain gates, eta_0, eta_0 / m, top-eigvec
    participation ratio, mean |Cor_ij| off the diagonal)."""
    S = np.where((p > lo) & (p < 1 - lo))[0]
    m = len(S)
    if m < 3:
        return m, np.nan, np.nan, np.nan, np.nan
    q = p[S]
    C = GG[np.ix_(S, S)] - np.outer(q, q)
    v = q * (1 - q)
    np.fill_diagonal(C, v)
    R_ = C / np.sqrt(np.outer(v, v))
    ev, evec = np.linalg.eigh(R_)
    u = evec[:, -1]
    off = np.abs(R_[~np.eye(m, dtype=bool)]).mean()
    return m, ev[-1] - 1.0, (ev[-1] - 1.0) / m, 1.0 / np.sum(u ** 4), off


def first_k(eps, ks, thr):
    for k, e in zip(ks, eps):
        if e <= thr:
            return k
    return np.nan


def corr_eps(eA, eB, nA, nB):
    """noise-corrected relative error: sqrt(<eA, eB> / <nA, nB>); eA, eB errors of the two atlases, nA, nB the normalisers."""
    num = float(np.sum(eA * eB))
    den = float(np.sum(nA * nB))
    return np.sqrt(max(num, 0.0) / den) if den > 0 else np.nan


# ----------------------------------------------------------------------------------------------------------- atlas loading
def load_atlas(path, source="ad"):
    """source 'ad': K_s = all-distinct part of kappa3(a_s) (default, the task's definition); 'full': the whole kappa3(a_s),
    slices included (its partial trace is non-zero, see the digest's section on the harmonic / trace split)."""
    z = np.load(path)
    W = z["weights"].astype(np.float64)
    L, n, _ = W.shape
    Phi = z["gate_p"].astype(np.float64)
    pre_s, post_s = z["pre_s"], z["post_s"]
    post_M11 = z["post_M11"].astype(np.float64)
    pre_M11 = z["pre_M11"].astype(np.float64)
    K, Ca, Cz, D = [], [], [], []
    P3 = z["post_M3"]
    for l in range(L):
        mua = post_s[0, l]
        k3 = central3(P3[l], post_M11[l], mua)
        K.append(all_distinct(k3) if source == "ad" else k3)
        Ca.append(post_M11[l] - np.outer(mua, mua))
    del P3
    Z3 = z["pre_M3"]
    for l in range(L):
        muz = pre_s[0, l]
        D.append(d21(central3(Z3[l], pre_M11[l], muz)))
        Cz.append(pre_M11[l] - np.outer(muz, muz))
    del Z3
    GG = z["gate_GG"].astype(np.float64) if "gate_GG" in z.files else None
    return dict(W=W, Phi=Phi, K=K, Ca=Ca, Cz=Cz, D=D, N=int(z["n_samples"]), seed=int(z["sample_seed"]), GG=GG,
                muz=np.array(pre_s[0], dtype=np.float64), mua=np.array(post_s[0], dtype=np.float64))


# ----------------------------------------------------------------------------------------------------------- mode: atlas
def run_atlas(pA, pB=None, save=None, kgrid=KGRID, rng_seed=0, source="ad"):
    t0 = time.time()
    a = load_atlas(pA, source)
    b = load_atlas(pB, source) if pB else None
    W = a["W"]
    L, n, _ = W.shape
    if b is not None:
        assert np.allclose(W, b["W"]), "atlas B must be the same MLP"
        Phi = 0.5 * (a["Phi"] + b["Phi"])
    else:
        Phi = a["Phi"]
    kgrid = [k for k in kgrid if k <= n]
    kc = [k for k in KCOARSE if k <= n]
    print(f"# transfer_spectrum atlas: {pA}" + (f" | {pB}" if pB else "") + f"; width {n}, depth {L}, N = {a['N']}"
          + (f" + {b['N']}" if b else "") + f"; source K_s = {'all-distinct part of' if source == 'ad' else 'whole'} kappa3(a_s)"
          + f"; load {time.time() - t0:.0f} s", flush=True)
    A = [W[l].T.copy() for l in range(L)]

    # ---------------- (a) propagators and their spectra
    J, SV = {}, {}
    for s in range(L - 1):
        M = A[s + 1].copy()
        J[(s, s + 1)] = M
        for t in range(s + 2, L):
            M = A[t] @ (Phi[t - 1][:, None] * M)
            J[(s, t)] = M
    for key, M in J.items():
        SV[key] = np.linalg.svd(M)
    rA = np.array([rmom(A[l]) for l in range(L)])
    rD = np.array([np.mean(Phi[l] ** 4) / np.mean(Phi[l] ** 2) ** 2 for l in range(L)])
    print("\n## (a) spectra of the propagators J_{s->t} (a_s -> z_t)")
    print("per-layer r(A^T A) [Ginibre: 2]: " + " ".join(f"{x:.2f}" for x in rA))
    print("per-layer r(D^2) = E Phi^4 / (E Phi^2)^2 [all-equal gates: 1; 0/1 mask with p = 1/2: 2]: " + " ".join(f"{x:.2f}" for x in rD))
    stats = {}
    for (s, t), (U, sv, Vt) in SV.items():
        pr, r90, r99, r999 = spec_stats(sv)
        r_pred = 1.0 + np.sum(rA[s + 1:t + 1] - 1.0) + np.sum(rD[s + 1:t] - 1.0)
        prG = spec_stats(np.linalg.svd(Phi[t][:, None] * J[(s, t)], compute_uv=False))[0]
        stats[(s, t)] = dict(pr=pr, r90=r90, r99=r99, r999=r999, pr_free=n / r_pred, prG=prG)
    print("\nby age a = t - s (mean over the L - a pairs): PR measured | PR_free = n / r_pred | r90 r99 r999 | PR of G = D_t J"
          " | s1/s2, s1/s8, s1/s32")
    for age in range(1, L):
        ks = [(s, s + age) for s in range(L - age)]
        m = lambda f: np.mean([stats[k][f] for k in ks])
        rat = np.mean([[SV[k][1][0] / SV[k][1][1], SV[k][1][0] / SV[k][1][7], SV[k][1][0] / SV[k][1][31]] for k in ks], 0)
        print(f"age {age:>2} ({len(ks):>2}): PR {m('pr'):6.1f} | free {m('pr_free'):6.1f} | r90 {m('r90'):5.1f} r99 {m('r99'):5.1f}"
              f" r999 {m('r999'):5.1f} | PR_G {m('prG'):6.1f} | {rat[0]:5.2f} {rat[1]:6.2f} {rat[2]:7.2f}")
    print("\nmean-direction outlier by age (mean over pairs): s1^2 / mean s^2 | s1/s2 | <u_1(J_{s->t}), mu_z(t)/|mu_z(t)|>^2 |"
          " ||J mu_a(s)||^2 / (|mu_a(s)|^2 mean s^2)")
    for age in range(1, L):
        ks = [(s, s + age) for s in range(L - age)]
        o = np.mean([[SV[k][1][0] ** 2 / np.mean(SV[k][1] ** 2), SV[k][1][0] / SV[k][1][1],
                      (SV[k][0][:, 0] @ a["muz"][k[1]]) ** 2 / np.sum(a["muz"][k[1]] ** 2),
                      np.sum((J[k] @ a["mua"][k[0]]) ** 2) / np.sum(a["mua"][k[0]] ** 2) / np.mean(SV[k][1] ** 2)] for k in ks], 0)
        print(f"age {age:>2}: {o[0]:6.1f} | {o[1]:4.2f} | {o[2]:4.2f} | {o[3]:5.2f}")
    print("\nPR matrix (rows s, columns t):")
    print("     " + "".join(f"{t:>6}" for t in range(1, L)))
    for s in range(L - 1):
        print(f"s={s:>2} " + "".join(f"{stats[(s, t)]['pr']:6.1f}" if t > s else "      " for t in range(1, L)))
    sv0 = SV[(0, L - 1)][1]
    with np.errstate(divide="ignore"):
        lyap = np.log(sv0) / (L - 1)
    print(f"\nfinite-time Lyapunov exponents of J_(0->{L - 1}) (log s_i / {L - 1}): top 8 " + " ".join(f"{x:+.3f}" for x in lyap[:8])
          + f" | i = 16, 32, 64, 96, 128: " + " ".join(f"{lyap[i - 1]:+.3f}" for i in (16, 32, 64, 96, n)))
    gaps = -np.diff(lyap)
    print(f"mean consecutive Lyapunov gap over i = 1..16: {gaps[:16].mean():.4f}; i = 17..64: {gaps[16:64].mean():.4f}"
          f" (an O(1) gap would be the 'expander' regime; 1/(2n) = {1 / (2 * n):.4f})")
    print("\nalignment of top-k singular subspaces, ||Q_k^T Q'_k||_F^2 / k (random: k/n):")
    for k in (8, 16, 32):
        for t in (L - 1, L // 2 + 2):
            row = []
            for s in range(t - 1):
                Uk, Uref = SV[(s, t)][0][:, :k], SV[(0, t)][0][:, :k]
                row.append(np.sum((Uk.T @ Uref) ** 2) / k)
            print(f"  U side (target t = {t:>2}, k = {k:>2}), source s = 0..{t - 2} vs s = 0: " + " ".join(f"{x:.2f}" for x in row))
        for s in (1, 5):
            row = []
            for t in range(s + 1, L):
                Vk, Vref = SV[(s, t)][2][:k].T, SV[(s, L - 1)][2][:k].T
                row.append(np.sum((Vk.T @ Vref) ** 2) / k)
            print(f"  V side (source s = {s}, k = {k:>2}), target t = {s + 1}..{L - 1} vs t = {L - 1}: " + " ".join(f"{x:.2f}" for x in row))
    print(f"  (random baseline k/n: " + ", ".join(f"k={k}: {k / n:.2f}" for k in (8, 16, 32)) + ")", flush=True)

    # ---------------- (b) the transported all-distinct third cumulant
    rng = np.random.default_rng(rng_seed)
    nk = len(kgrid)
    keys = [(s, t) for s in range(L - 1) for t in range(s + 1, L)]
    idx = {k: i for i, k in enumerate(keys)}
    P = len(keys)
    R = dict(
        share=np.full(P, np.nan), share_raw=np.full(P, np.nan), unexpl=np.full(P, np.nan), unexpl_raw=np.full(P, np.nan),
        own=np.full((P, nk), np.nan), tot=np.full((P, nk), np.nan),            # noise-corrected (or raw if no B)
        own_raw=np.full((P, nk), np.nan), tot_raw=np.full((P, nk), np.nan),
        own_ens=np.full((P, nk), np.nan), noise_d=np.full(P, np.nan))
    BASES = ("rand", "hosvd", "fix", "covz", "covprop")
    for nm in BASES:
        R["own_" + nm] = np.full((P, len(kc)), np.nan)
        R["tot_" + nm] = np.full((P, len(kc)), np.nan)
    Q = {}       # source-side fixed basis V(J_{s->L-1}) and the source cores in it
    for s in range(L - 1):
        V = SV[(s, L - 1)][2].T
        Q[s] = (V, rot3(a["K"][s], V), rot3(b["K"][s], V) if b else None)
    cur = {}
    print("\n## (b) all-distinct kappa3(a_s) transported by J_{s->t}^(x)3 (births ignored); D21 at z_t")
    for t in range(1, L):
        g3 = None
        nxt = {}
        for s, (XA, XB) in cur.items():
            g = Phi[t - 1]
            nxt[s] = (rot3(phi3(XA, g), W[t]), rot3(phi3(XB, g), W[t]) if XB is not None else None)
        nxt[t - 1] = (rot3(a["K"][t - 1], W[t]), rot3(b["K"][t - 1], W[t]) if b else None)
        cur = nxt
        DtA = a["D"][t]
        DtB = b["D"][t] if b else None
        for s in sorted(cur):
            i = idx[(s, t)]
            XA, XB = cur[s]
            U, sv, Vt = SV[(s, t)]
            dA = d21(XA)
            dB = d21(XB) if XB is not None else None
            R["share_raw"][i] = np.linalg.norm(dA) / np.linalg.norm(DtA)
            R["unexpl_raw"][i] = np.linalg.norm(DtA - dA) / np.linalg.norm(DtA)
            if b is not None:
                R["share"][i] = np.sqrt(max(np.sum(dA * dB), 0) / np.sum(DtA * DtB))
                R["unexpl"][i] = corr_eps(DtA - dA, DtB - dB, DtA, DtB)
                R["noise_d"][i] = np.linalg.norm(dA - dB) / np.sqrt(2) / np.sqrt(max(np.sum(dA * dB), 1e-300))
            YA = rot3(XA, U)
            YB = rot3(XB, U) if XB is not None else None
            for j, k in enumerate(kgrid):
                eA = d21_from_core(YA[:k, :k, :k], U[:, :k]) - dA
                R["own_raw"][i, j] = np.linalg.norm(eA) / np.linalg.norm(dA)
                R["tot_raw"][i, j] = np.linalg.norm(eA) / np.linalg.norm(DtA)
                if YB is not None:
                    eB = d21_from_core(YB[:k, :k, :k], U[:, :k]) - dB
                    R["own"][i, j] = corr_eps(eA, eB, dA, dB)
                    R["tot"][i, j] = corr_eps(eA, eB, DtA, DtB)
            if b is None:
                R["own"][i], R["tot"][i] = R["own_raw"][i], R["tot_raw"][i]
            ens = eps_ensemble(U, sv)
            R["own_ens"][i] = ens[np.array(kgrid) - 1]
            # baselines (atlas A, own-relative; noise-corrected when B is present)
            Rq = np.linalg.qr(rng.standard_normal((n, n)))[0]
            G = XA.reshape(n, n * n)
            ev, H = np.linalg.eigh(G @ G.T)
            H = H[:, ::-1]
            YR, YH = rot3(XA, Rq), rot3(XA, H)
            YRb = rot3(XB, Rq) if XB is not None else None
            YHb = rot3(XB, H) if XB is not None else None
            V, CA, CB = Q[s]
            JV = J[(s, t)] @ V
            # bases a chain already has: eigenvectors of the current covariance C(z_t), and of the transported covariance
            Cz = a["Cz"][t]
            Ez = np.linalg.eigh(0.5 * (Cz + Cz.T))[1][:, ::-1]
            cp = J[(s, t)] @ a["Ca"][s] @ J[(s, t)].T
            Ep = np.linalg.eigh(0.5 * (cp + cp.T))[1][:, ::-1]
            YZ, YP = rot3(XA, Ez), rot3(XA, Ep)
            YZb = rot3(XB, Ez) if XB is not None else None
            YPb = rot3(XB, Ep) if XB is not None else None
            for j, k in enumerate(kc):
                for name, YY, YYb, B_ in (("rand", YR, YRb, Rq[:, :k]), ("hosvd", YH, YHb, H[:, :k]),
                                          ("fix", CA, CB, JV[:, :k]), ("covz", YZ, YZb, Ez[:, :k]), ("covprop", YP, YPb, Ep[:, :k])):
                    eA = d21_from_core(YY[:k, :k, :k], B_) - dA
                    if YYb is not None:
                        eB = d21_from_core(YYb[:k, :k, :k], B_) - dB
                        R["own_" + name][i, j] = corr_eps(eA, eB, dA, dB)
                        R["tot_" + name][i, j] = corr_eps(eA, eB, DtA, DtB)
                    else:
                        R["own_" + name][i, j] = np.linalg.norm(eA) / np.linalg.norm(dA)
                        R["tot_" + name][i, j] = np.linalg.norm(eA) / np.linalg.norm(DtA)
        print(f"  t = {t:>2} done ({time.time() - t0:.0f} s)", flush=True)

    tag = "noise-corrected with the second atlas" if b else "raw (single atlas)"
    print(f"\nper pair (s, t): share = ||d|| / ||D21(z_t)||, unexpl = ||D21(z_t) - d|| / ||D21(z_t)||, then eps_tot(k) [{tag}]")
    show = [4, 8, 16, 32, 64]
    jj = [kgrid.index(k) for k in show if k in kgrid]
    for (s, t) in keys:
        i = idx[(s, t)]
        print(f"s={s:>2} t={t:>2} age={t - s:>2} share {R['share'][i] if b else R['share_raw'][i]:5.3f} unexpl "
              f"{R['unexpl'][i] if b else R['unexpl_raw'][i]:5.3f} | eps_own " + " ".join(f"k{kgrid[j]}:{R['own'][i, j]:.3f}" for j in jj)
              + " | eps_tot " + " ".join(f"k{kgrid[j]}:{R['tot'][i, j]:.3f}" for j in jj))
    print(f"\nby age (median over pairs; [{tag}]): share, unexpl | eps_own(k) | eps_tot(k) | k_2% own (median, max) | k_2% tot (median, max)")
    for age in range(1, L):
        ii = [idx[(s, s + age)] for s in range(L - age)]
        med = lambda x: np.nanmedian(x[ii], 0)
        k2o = [first_k(R["own"][i], kgrid, 0.02) for i in ii]
        k2t = [first_k(R["tot"][i], kgrid, 0.02) for i in ii]
        sh = R["share"] if b else R["share_raw"]
        ue = R["unexpl"] if b else R["unexpl_raw"]
        print(f"age {age:>2}: share {np.nanmedian(sh[ii]):.3f} unexpl {np.nanmedian(ue[ii]):.3f} | own "
              + " ".join(f"{med(R['own'])[j]:.3f}" for j in jj) + " | tot " + " ".join(f"{med(R['tot'])[j]:.3f}" for j in jj)
              + f" | k2own {np.nanmedian(k2o):5.1f} {np.nanmax(k2o):5.0f} | k2tot {np.nanmedian(k2t):5.1f} {np.nanmax(k2t):5.0f}")
    print(f"\nbaselines by age (median eps_own over pairs) at k = " + ", ".join(map(str, kc)))
    jt = [kgrid.index(k) for k in kc]
    for age in range(1, L):
        ii = [idx[(s, s + age)] for s in range(L - age)]
        med = lambda x: np.nanmedian(x[ii], 0)
        print(f"age {age:>2}: top-k of J " + " ".join(f"{x:.3f}" for x in med(R['own'])[jt]) + " | random " + " ".join(f"{x:.3f}" for x in med(R['own_rand']))
              + " | HOSVD of X " + " ".join(f"{x:.3f}" for x in med(R['own_hosvd'])) + " | fixed V_k(J_{s->L-1}) at birth "
              + " ".join(f"{x:.3f}" for x in med(R['own_fix'])) + " | eig C(z_t) " + " ".join(f"{x:.3f}" for x in med(R['own_covz']))
              + " | eig J C(a_s) J^T " + " ".join(f"{x:.3f}" for x in med(R['own_covprop'])))
    print("\nGaussian-ensemble prediction (d) vs measured eps_own (median over pairs) at k = " + ", ".join(map(str, show)))
    for age in range(1, L):
        ii = [idx[(s, s + age)] for s in range(L - age)]
        med = lambda x: np.nanmedian(x[ii], 0)
        print(f"age {age:>2}: measured " + " ".join(f"{med(R['own'])[j]:.3f}" for j in jj) + " | ensemble "
              + " ".join(f"{med(R['own_ens'])[j]:.3f}" for j in jj)
              + f" | k2own measured {np.nanmedian([first_k(R['own'][i], kgrid, 0.02) for i in ii]):5.1f}"
              + f" ensemble {np.nanmedian([first_k(R['own_ens'][i], kgrid, 0.02) for i in ii]):5.1f}")
    print("\nold tier as one tensor: target t, age threshold w -> X_{t-w,t} carries everything present at layer t-w; k_2% of eps_tot"
          " with the bases U_k(J_{t-w->t}) / birth-fixed V_k(J_{t-w->L-1}) / eig C(z_t) / HOSVD of X (the last three on the coarse grid)")
    print("      " + "".join(f"  w={w:<14}" for w in range(2, 9)))
    for t in range(3, L):
        row = []
        for w in range(2, 9):
            if t - w < 0:
                row.append("       -        ")
                continue
            i = idx[(t - w, t)]
            ks_ = [first_k(R["tot"][i], kgrid, 0.02)] + [first_k(R["tot_" + nm][i], kc, 0.02) for nm in ("fix", "covz", "hosvd")]
            row.append("/".join(f"{x:.0f}" for x in ks_).ljust(16))
        print(f"t={t:>2} " + " ".join(row))
    print("\nk_2% of eps_own by age (median over pairs), per basis: U_k(J) | fixed | eig C(z_t) | eig J C(a_s) J^T | HOSVD")
    for age in range(1, L):
        ii = [idx[(s, s + age)] for s in range(L - age)]
        vals = [np.nanmedian([first_k(R["own"][i], kgrid, 0.02) for i in ii])]
        vals += [np.nanmedian([first_k(R["own_" + nm][i], kc, 0.02) for i in ii]) for nm in ("fix", "covz", "covprop", "hosvd")]
        print(f"age {age:>2}: " + " | ".join(f"{v:5.1f}" for v in vals))

    # ---------------- (c) the symmetric-matrix transfer map
    print("\n## (c) S -> B S B^T on Sym_n, B_l = A_{l+1} D_l (z_l -> z_{l+1}); singular values s_i s_j, eigenvalues l_i l_j (i <= j)")
    print("layer: op norm s1^2 | spectral radius |l1|^2 | PR(Sym^2) | r90 r99 of the Sym^2 singular values | dim Sym^2 = n(n+1)/2 = "
          f"{n * (n + 1) // 2}")
    iu = np.triu_indices(n)
    def sym2_stats(sv):
        p = (sv[:, None] * sv[None, :])[iu]
        return spec_stats(np.sort(p)[::-1])
    for l in range(L - 1):
        B_ = A[l + 1] * Phi[l][None, :]
        sv = np.linalg.svd(B_, compute_uv=False)
        lam = np.abs(np.linalg.eigvals(B_))
        pr, r90, r99, _ = sym2_stats(sv)
        print(f"  l={l:>2}: {sv[0] ** 2:7.3f} | {lam.max() ** 2:6.3f} | {pr:8.1f} | {r90:5d} {r99:5d}")
    print("products z_s -> z_t (B_{t-1} ... B_s = J_{s->t} D_s), by age: PR(Sym^2), r90, r99 (mean over pairs)")
    for age in range(1, L):
        rows = []
        for s in range(L - age):
            sv = np.linalg.svd(J[(s, s + age)] * Phi[s][None, :], compute_uv=False)
            rows.append(sym2_stats(sv)[:3])
        r = np.mean(rows, 0)
        print(f"  age {age:>2}: PR {r[0]:8.1f} r90 {r[1]:6.1f} r99 {r[2]:6.1f}")
    print("transported covariance c = J_{s->t} C(a_s) J^T: share ||c|| / ||C(z_t)||, and ||P_k c P_k - c|| / ||c|| with P_k = U_k U_k^T, by age (median)")
    for age in range(1, L):
        sh, ek = [], []
        for s in range(L - age):
            t = s + age
            Jm = J[(s, t)]
            c = Jm @ a["Ca"][s] @ Jm.T
            U = SV[(s, t)][0]
            sh.append(np.linalg.norm(c) / np.linalg.norm(a["Cz"][t]))
            ek.append([np.linalg.norm(U[:, :k] @ (U[:, :k].T @ c @ U[:, :k]) @ U[:, :k].T - c) / np.linalg.norm(c) for k in show])
        print(f"  age {age:>2}: share {np.median(sh):.3f} | " + " ".join(f"k{k}:{x:.4f}" for k, x in zip(show, np.median(ek, 0))))
    if a["GG"] is not None:
        print("\n## (e) the gate (face) law of each layer as a Markov random field on {0,1}^n: unpinned spectral independence")
        print("layer: m uncertain gates (0.02 < p < 0.98) | eta_0 = lambda_max(Psi) | eta_0 / m | PR of the top eigenvector | mean |Cor_ij|"
              " | same for the Gaussian surrogate (gates of N(mu, C(z_l)) via Sheppard's arcsine law)")
        for l in range(L):
            p_ = a["Phi"][l]
            m_, e0, e0m, prv, off = gate_si(a["GG"][l], p_)
            # Gaussian surrogate: P(z_i > 0, z_j > 0) for a bivariate normal is not closed-form with means; use the zero-mean
            # (Sheppard) form on standardised correlations as a reference only
            Cz = a["Cz"][l]
            sd = np.sqrt(np.diag(Cz))
            rho = np.clip(Cz / np.outer(sd, sd), -1, 1)
            GGs = 0.25 + np.arcsin(rho) / (2 * np.pi)
            ms, es, esm, prs, offs = gate_si(GGs, np.full(n, 0.5))
            print(f"  l={l:>2}: {m_:4d} | {e0:7.2f} | {e0m:.3f} | {prv:6.1f} | {off:.3f} || Sheppard: {es:7.2f} | {esm:.3f} | {prs:6.1f} | {offs:.3f}")
    print(f"\n# done in {time.time() - t0:.0f} s", flush=True)
    if save:
        out = {k: v for k, v in R.items()}
        out.update(keys=np.array(keys), kgrid=np.array(kgrid), kcoarse=np.array(kc), rA=rA, rD=rD,
                   pr=np.array([stats[k]["pr"] for k in keys]), pr_free=np.array([stats[k]["pr_free"] for k in keys]),
                   r90=np.array([stats[k]["r90"] for k in keys]), r99=np.array([stats[k]["r99"] for k in keys]),
                   r999=np.array([stats[k]["r999"] for k in keys]), sv=np.array([SV[k][1] for k in keys]))
        np.savez_compressed(save, **out)


# ----------------------------------------------------------------------------------------------------------- mode: scaling
def mc_gates(W, N, rng, chunk=4000):
    """gate probabilities P(z_l > 0) and the means of z_l and a_l from N Monte Carlo inputs x ~ N(0, I)."""
    L, n, _ = W.shape
    cnt, sz, sa = np.zeros((L, n)), np.zeros((L, n)), np.zeros((L, n))
    GG = np.zeros((L, n, n))
    done = 0
    while done < N:
        m = min(chunk, N - done)
        h = rng.standard_normal((m, n))
        for l in range(L):
            z = h @ W[l]
            g = (z > 0).astype(np.float64)
            cnt[l] += g.sum(0)
            GG[l] += g.T @ g
            sz[l] += z.sum(0)
            h = np.maximum(z, 0.0)
            sa[l] += h.sum(0)
        done += m
    return cnt / N, sz / N, sa / N, GG / N


def run_scaling(widths, L, N, seed, save=None):
    t0 = time.time()
    fr = np.array([1 / 64, 1 / 32, 1 / 16, 1 / 8, 3 / 16, 1 / 4, 3 / 8, 1 / 2, 5 / 8, 3 / 4, 7 / 8, 1.0])
    print(f"# transfer_spectrum scaling: random He MLPs (W ~ N(0, 2/n), no bias, x ~ N(0, I)), depth {L}, gates from {N} Monte Carlo"
          f" samples; Gaussian-ensemble eps_own(k) (part (d))", flush=True)
    res = {}
    for n in widths:
        rng = np.random.default_rng(seed + n)
        W = rng.standard_normal((L, n, n)) * np.sqrt(2.0 / n)
        Phi, muz, mua, GGm = mc_gates(W, N, rng)
        si = [gate_si(GGm[l], Phi[l]) for l in range(L)]
        rD = np.array([np.mean(Phi[l] ** 4) / np.mean(Phi[l] ** 2) ** 2 for l in range(L)])
        A = [W[l].T.copy() for l in range(L)]
        rA = np.array([rmom(A[l]) for l in range(L)])
        rows = {}
        for s in range(L - 1):
            M = A[s + 1].copy()
            for t in range(s + 1, L):
                if t > s + 1:
                    M = A[t] @ (Phi[t - 1][:, None] * M)
                U, sv, _ = np.linalg.svd(M)
                pr, r90, r99, r999 = spec_stats(sv)
                r_pred = 1.0 + np.sum(rA[s + 1:t + 1] - 1.0) + np.sum(rD[s + 1:t] - 1.0)
                ens = eps_ensemble(U, sv)
                ks = np.arange(1, n + 1)
                mz = muz[t] / np.linalg.norm(muz[t])
                ma = mua[s] / np.linalg.norm(mua[s])
                rows[(s, t)] = dict(pr=pr, r90=r90, r99=r99, r999=r999, pr_free=n / r_pred,
                                    out=sv[0] ** 2 / np.mean(sv ** 2), s12=sv[0] / sv[1], u1mu=float((U[:, 0] @ mz) ** 2),
                                    gmu=float(np.linalg.norm(M @ ma) ** 2 / np.mean(sv ** 2)),
                                    k1=first_k(ens, ks, 0.01), k2=first_k(ens, ks, 0.02), k5=first_k(ens, ks, 0.05),
                                    k10=first_k(ens, ks, 0.10), k20=first_k(ens, ks, 0.20),
                                    efr=ens[np.maximum(np.round(fr * n).astype(int), 1) - 1])
            print(f"  n = {n}: source s = {s} done ({time.time() - t0:.0f} s)", flush=True)
        res[n] = dict(rows=rows, rD=rD, rA=rA, si=si)
        print(f"\nn = {n}: r(D^2) per layer " + " ".join(f"{x:.2f}" for x in rD))
        print(f"n = {n}, by age (mean over pairs): PR | PR_free | PR/n | r90/n r99/n r999/n | ensemble k for eps_own = 20, 10, 5, 2, 1 %"
              f" (as k and as k/n)")
        for age in range(1, L):
            ks_ = [(s, s + age) for s in range(L - age)]
            m = lambda f: np.mean([rows[k][f] for k in ks_])
            print(f"  age {age:>2}: PR {m('pr'):7.1f} | free {m('pr_free'):7.1f} | {m('pr') / n:.3f} | {m('r90') / n:.3f} {m('r99') / n:.3f}"
                  f" {m('r999') / n:.3f} | k {m('k20'):6.1f} {m('k10'):6.1f} {m('k5'):6.1f} {m('k2'):6.1f} {m('k1'):6.1f}"
                  f" | k/n {m('k20') / n:.3f} {m('k10') / n:.3f} {m('k5') / n:.3f} {m('k2') / n:.3f} {m('k1') / n:.3f}", flush=True)
    print("\ngate (face) law, unpinned spectral independence by layer and width: m uncertain gates, eta_0, eta_0 / m"
          " (eta_0 = lambda_max(Psi), Monte Carlo pair gate statistics)")
    for l in range(L):
        print(f"l={l:>2}  " + "  ".join(f"n={n}: {res[n]['si'][l][0]:4d} {res[n]['si'][l][1]:7.2f} {res[n]['si'][l][2]:.3f}" for n in widths))
    print("\nsummary: ensemble k_2% / n by age and width (mean over pairs)")
    print("age  " + "".join(f"{n:>9}" for n in widths))
    for age in range(1, L):
        print(f"{age:>3}  " + "".join(f"{np.mean([res[n]['rows'][(s, s + age)]['k2'] for s in range(L - age)]) / n:9.3f}" for n in widths))
    print("\nmean-direction outlier, by age and width (mean over pairs): s1^2 / mean s^2 | s1/s2 | <u_1(J_{s->t}), mu_z(t)/|mu_z(t)|>^2 |"
          " ||J mu_a(s)||^2 / (|mu_a(s)|^2 mean s^2)")
    for age in range(1, L):
        print(f"{age:>3}  " + "  ".join(
            f"n={n}: {np.mean([res[n]['rows'][(s, s + age)]['out'] for s in range(L - age)]):6.1f} "
            f"{np.mean([res[n]['rows'][(s, s + age)]['s12'] for s in range(L - age)]):4.2f} "
            f"{np.mean([res[n]['rows'][(s, s + age)]['u1mu'] for s in range(L - age)]):4.2f} "
            f"{np.mean([res[n]['rows'][(s, s + age)]['gmu'] for s in range(L - age)]):5.1f}" for n in widths))
    print("\nsummary: PR * age / n by age and width (free probability predicts ~ age / (2 age - 1 + ...) for gated products)")
    print("age  " + "".join(f"{n:>9}" for n in widths))
    for age in range(1, L):
        print(f"{age:>3}  " + "".join(f"{np.mean([res[n]['rows'][(s, s + age)]['pr'] for s in range(L - age)]) * age / n:9.3f}" for n in widths))
    print(f"\n# done in {time.time() - t0:.0f} s", flush=True)
    if save:
        out = {}
        for n in widths:
            keys = sorted(res[n]["rows"])
            for f in ("pr", "pr_free", "r90", "r99", "r999", "k1", "k2", "k5", "k10", "k20", "out", "s12", "u1mu", "gmu"):
                out[f"n{n}_{f}"] = np.array([res[n]["rows"][k][f] for k in keys])
            out[f"n{n}_efr"] = np.array([res[n]["rows"][k]["efr"] for k in keys])
            out[f"n{n}_keys"] = np.array(keys)
            out[f"n{n}_rD"] = res[n]["rD"]
        out["fractions"] = fr
        np.savez_compressed(save, **out)


# ----------------------------------------------------------------------------------------------------------- mode: selftest
def selftest():
    rng = np.random.default_rng(1)
    ok = True
    # rot3 and d21_from_core
    n, m = 7, 5
    X = rng.standard_normal((n, n, n))
    U = rng.standard_normal((n, m))
    e1 = np.max(np.abs(rot3(X, U) - np.einsum("ijk,ip,jq,kr->pqr", X, U, U, U)))
    Y = rng.standard_normal((m, m, m))
    e2 = np.max(np.abs(d21_from_core(Y, U) - d21(np.einsum("pqr,ap,bq,cr->abc", Y, U, U, U))))
    print(f"rot3 vs einsum {e1:.1e}; d21_from_core vs explicit {e2:.1e}")
    ok &= e1 < 1e-10 and e2 < 1e-10
    # Sym^2 functor: singular values and eigenvalues of S -> B S B^T on Sym_n
    n = 9
    B = rng.standard_normal((n, n)) * (rng.random(n) > 0.4)[None, :]
    basis = []
    for i in range(n):
        for j in range(i, n):
            E = np.zeros((n, n))
            if i == j:
                E[i, i] = 1.0
            else:
                E[i, j] = E[j, i] = 1 / np.sqrt(2)
            basis.append(E)
    T = np.array([[np.sum(Eo * (B @ Ei @ B.T)) for Ei in basis] for Eo in basis])
    sv = np.linalg.svd(B, compute_uv=False)
    lam = np.linalg.eigvals(B)
    iu = np.triu_indices(n)
    p_sv = np.sort((sv[:, None] * sv[None, :])[iu])[::-1]
    p_ev = np.sort(np.abs((lam[:, None] * lam[None, :])[iu]))[::-1]
    e3 = np.max(np.abs(np.linalg.svd(T, compute_uv=False) - p_sv))
    e4 = np.max(np.abs(np.sort(np.abs(np.linalg.eigvals(T)))[::-1] - p_ev))
    print(f"Sym^2: singular values = s_i s_j (i <= j) to {e3:.1e}; |eigenvalues| = |l_i l_j| to {e4:.1e}")
    ok &= e3 < 1e-8 and e4 < 1e-6
    # Gaussian-ensemble formula against Monte Carlo (symmetric Gaussian tensors, and their all-distinct parts)
    n, reps = 32, 300
    Jm = rng.standard_normal((n, n)) * np.sqrt(2.0 / n)
    for _ in range(2):
        Jm = (rng.standard_normal((n, n)) * np.sqrt(2.0 / n)) @ ((rng.random(n) < 0.7)[:, None] * Jm)
    U, sv, _ = np.linalg.svd(Jm)
    ens = eps_ensemble(U, sv)
    for k in (2, 4, 8, 12, 16):
        num = {"sym": 0.0, "ad": 0.0}
        den = {"sym": 0.0, "ad": 0.0}
        for r in range(reps):
            Z = rng.standard_normal((n, n, n))
            Ks = (Z + Z.transpose(0, 2, 1) + Z.transpose(1, 0, 2) + Z.transpose(1, 2, 0) + Z.transpose(2, 0, 1) + Z.transpose(2, 1, 0)) / 6
            for name, K in (("sym", Ks), ("ad", all_distinct(Ks))):
                Xt = rot3(K, Jm.T)
                d = d21(Xt)
                dk = d21_from_core(rot3(Xt, U[:, :k]), U[:, :k])
                num[name] += np.sum((dk - d) ** 2)
                den[name] += np.sum(d ** 2)
        print(f"ensemble formula k = {k:>2}: predicted {ens[k - 1]:.3f} | Monte Carlo symmetric {np.sqrt(num['sym'] / den['sym']):.3f}"
              f" | all-distinct {np.sqrt(num['ad'] / den['ad']):.3f}")
        ok &= abs(ens[k - 1] - np.sqrt(num["sym"] / den["sym"])) < 0.03
    # free-probability PR of masked Gaussian products
    n = 400
    Jm = np.eye(n)
    r_pred = 1.0
    for l in range(6):
        p = 0.5
        mask = (rng.random(n) < p).astype(float)
        Al = rng.standard_normal((n, n)) / np.sqrt(n)
        Jm = Al @ (mask[:, None] * Jm) if l else Al
        r_pred += (rmom(Al) - 1.0) + ((1.0 / p - 1.0) if l else 0.0)
    pr = spec_stats(np.linalg.svd(Jm, compute_uv=False))[0]
    print(f"free probability, product of 6 Gaussian 400 x 400 with 5 Bernoulli(1/2) masks: PR {pr:.1f} vs n / r_pred {n / r_pred:.1f}"
          f" (Ginibre-only would be n / 7 = {n / 7:.1f})")
    ok &= abs(pr - n / r_pred) / (n / r_pred) < 0.15
    print("SELFTEST", "PASSED" if ok else "FAILED")
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mode", choices=["atlas", "scaling", "selftest"])
    ap.add_argument("paths", nargs="*")
    ap.add_argument("--save", default=None)
    ap.add_argument("--widths", default="128,256,512,1024")
    ap.add_argument("--depth", type=int, default=16)
    ap.add_argument("--mc", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--source", choices=["ad", "full"], default="ad")
    args = ap.parse_args()
    if args.mode == "selftest":
        sys.exit(0 if selftest() else 1)
    if args.mode == "atlas":
        run_atlas(args.paths[0], args.paths[1] if len(args.paths) > 1 else None, save=args.save, source=args.source)
    else:
        run_scaling([int(x) for x in args.widths.split(",")], args.depth, args.mc, args.seed, save=args.save)


if __name__ == "__main__":
    main()
