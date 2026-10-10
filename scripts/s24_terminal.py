"""Terminal means: the three decisive measurements for a second pass that responds to the first pass's final error.

  python scripts/s24_terminal.py DATA OUT NET [NET ...]        (V56=dir with nXXX/out_X.npy, optional)

Estimator under test (terminal-means note). Angular frame: F is degree-1 homogeneous, x = R theta, R ~ chi_d,
theta ~ sigma uniform, so mu = E F(x) = int Y dsigma with Y(theta) = E[R] F(theta). With a first-pass answer b,
an orthonormal response basis V (n x r, P = V V^T), a control c of known integral a and a sampling density q,
    mu_hat = b + V [a - V^T b + (1/M) sum_j z(theta_j) / q(theta_j)],      z = V^T Y - c,
    risk   = |(I - P)(mu - b)|^2 / n  +  (J(q) - |zbar|^2) / (n M),        J(q) = int |z|^2 / q dsigma.
At the 0.1 B floor (204.8 n^3) after a CC1-sized pass 1 (~55 n^3) one sample costs 32.8 n^2, so M = 4500.
The floor target is risk < 3.1e-8 (adjusted 3.1e-9 = v56).

(1) best correction inside a response space: captured fraction |P e|^2/|e|^2, e = mu - b, for b = Gaussian closure,
    CC1, v56; families: quiet eigenspaces of the output covariance (pilot and pass-1), coherent functions of
    a = W_L mu_hat, coherent functions of every layer carried to the output by the closure tangent, left singular
    subspaces of the closure tangent T_(L<-l), final-layer readout-shape features, and the CC1 - closure direction.
(2) actual projected residual covariance on held-out samples, raw and with exact-mean controls fitted on the pilot
    (theta; theta and h_1 - E h_1; and the first-layer Gram form theta^T (W_1^T W_1) theta - tr/d).
(3) held-out int |z|^2 Psi_H dsigma for metric codes H (H^2 = I, tr H = 0, H = 2P_S - I): S = top half of
    W_1^T W_1, of the input-Jacobian Gram of the response space, of the pilot M_0; and the exact held-out
    IS second moment J(q_t) of the two-ellipsoid mixture and of the single ellipsoid along the pilot M_0.
Pilot set A and held-out set B are independent angular samples (2^16 each)."""
import sys, os, time, math, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from whest import gate_channel as gc
from whest.cc1 import CC1

BS, NA, NB, M_FLOOR, TARGET = 4096, 65536, 65536, 4500, 3.1e-8
rms = lambda x: float(np.sqrt(np.mean(np.square(x))))


def orth(X, tol=1e-10):
    """Orthonormal basis of the column span of X (rank-revealing via SVD)."""
    U, s, _ = np.linalg.svd(X, full_matrices=False)
    return U[:, s > tol * s[0]]


def mc(W32, N, seed):
    """Angular samples: theta (N, d) and Y = E[R] F(theta) (N, n), float32."""
    d = W32[0].shape[1]; n = W32[-1].shape[0]
    ER = math.sqrt(2.0) * math.exp(math.lgamma((d + 1) / 2) - math.lgamma(d / 2))
    Y = np.empty((N, n), np.float32); TH = np.empty((N, d), np.float32)
    for b in range(N // BS):
        g = np.random.default_rng([seed, b]).standard_normal((d, BS), dtype=np.float32)
        th = g / np.linalg.norm(g, axis=0, keepdims=True); h = th
        for Wl in W32: h = np.maximum(Wl @ h, 0)
        Y[b * BS:(b + 1) * BS] = (h * np.float32(ER)).T; TH[b * BS:(b + 1) * BS] = th.T
    return TH, Y, ER


def cov(Y, mean=None):
    m = Y.mean(0, dtype=np.float64) if mean is None else mean
    S = (Y.T @ Y).astype(np.float64) / len(Y)
    return S - np.outer(m, m), m


def hermite_feats(a, J):
    x = (a - a.mean()) / a.std(); H = [np.ones_like(x), x]
    for j in range(2, J): H.append(x * H[-1] - (j - 1) * H[-2])
    return np.stack(H[:J], 1)


def pass1(W64):
    t0 = time.time(); st = gc.closure_states(W64, keep_C=False)
    cl = np.array([s["m"] for s in st]); t1 = time.time()
    c = CC1(Q=7); cc, _ = c.run(W64); t2 = time.time()
    return st, cl, cc, c.last["K"], (t1 - t0, t2 - t1)


def tangents(W64, st):
    """Closure mean tangent T_(L<-l) = prod_(k>l) diag(Phi(alpha_k)) W_k, for l = -1 (input), 0, ..., L-1."""
    L = len(W64); T = {L - 1: np.eye(W64[-1].shape[0])}
    for l in range(L - 2, -2, -1):
        T[l] = (T[l + 1] * st[l + 1]["Pa"][None, :]) @ W64[l + 1]
    return T


def measure(D, OUT, net, V56):
    W32 = [np.ascontiguousarray(w, dtype=np.float32) for w in np.load(f"{D}/W_off{net}.npy")]
    W64 = [w.astype(np.float64) for w in W32]; L, n = len(W64), W64[0].shape[0]; d = W64[0].shape[1]
    T_ = np.load(f"{D}/truth_off{net}.npz"); truth = T_["m"].astype(np.float64); mu = truth[-1]
    st, cl, cc, Kcc, tp = pass1(W64)
    bs = {"closure": cl[-1], "CC1": cc[-1]}
    f56 = f"{V56}/n{net:03d}/out_{net}.npy" if V56 else ""
    if f56 and os.path.exists(f56): bs["v56"] = np.load(f56)[-1]
    E = {k: mu - v for k, v in bs.items()}
    t0 = time.time(); THA, YA, ER = mc(W32, NA, 24_000 + net); THB, YB, _ = mc(W32, NB, 25_000 + net); tmc = time.time() - t0
    SA, mA = cov(YA); SB, mB = cov(YB)
    lines = [f"\n=== net {net}  (pass 1: closure {tp[0]:.0f}s, CC1 {tp[1]:.0f}s; MC {NA + NB} angular samples {tmc:.0f}s) ===",
             f"first-pass final RMS: " + ", ".join(f"{k} {rms(e):.3e}" for k, e in E.items()) +
             f" | MC means vs truth: pilot {rms(mA - mu):.2e}, held-out {rms(mB - mu):.2e} (expected {math.sqrt(np.trace(SB) / n / NB):.2e})",
             f"output covariance: tr/n = {np.trace(SB) / n:.4f} (truth file avg_variance {float(T_['avg_variance']):.4f}; the angular frame removes the radial part)"]
    # ---------------- (1) best correction inside a response space ----------------
    lamA, UA = np.linalg.eigh(SA); UA = UA[:, ::-1]; lamA = lamA[::-1]
    lamB_on_A = np.sum(UA * (SB @ UA), 0)                             # held-out variance along pilot eigvecs
    pr = lambda lam: float(lam.sum() ** 2 / np.sum(lam ** 2))
    lamK, UK = np.linalg.eigh(Kcc); UK = UK[:, ::-1]; lamK = lamK[::-1]
    lamB_on_K = np.sum(UK * (SB @ UK), 0)
    lines.append(f"spectrum (pilot): PR {pr(lamA):.1f}; top-k share of trace k=1/4/16/64/256: " +
                 "/".join(f"{lamA[:k].sum() / lamA.sum():.3f}" for k in (1, 4, 16, 64, 256)) +
                 f"; bottom half carries {lamA[n // 2:].sum() / lamA.sum():.4f}; CC1 pass-1 K: PR {pr(np.clip(lamK, 0, None)):.1f}, "
                 f"tr K/tr S {lamK.sum() / np.trace(SB):.3f}")
    lines.append("(1a) quiet subspaces P_k = drop the top-k eigvecs; per b: captured fraction | risk at M=4500 "
                 "(= missed/n + held-out var/(nM)); basis from the pilot covariance / from CC1's pass-1 K")
    ks = [0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 768, 896, 960, 992, 1008, 1016, 1020]
    best = {}
    for name, U, lamB in (("pilot", UA, lamB_on_A), ("CC1-K", UK, lamB_on_K)):
        for bn, e in E.items():
            ec = (U.T @ e) ** 2; cume = np.concatenate([[0], np.cumsum(ec)]); tailv = np.concatenate([np.cumsum(lamB[::-1])[::-1], [0]])
            R = cume / n + tailv / (n * M_FLOOR); kb = int(np.argmin(R))
            best[(name, bn)] = (kb, R[kb], 1 - cume[kb] / cume[-1])
            orc = np.sum(np.minimum(ec, lamB / M_FLOOR)) / n
            lines.append(f"   {name:6s} b={bn:8s}: " + " ".join(f"k{k}:{1 - cume[k] / cume[-1]:.4f}|{R[k]:.1e}" for k in ks if k < n) +
                         f"\n{'':22s}best k={kb}: captured {1 - cume[kb] / cume[-1]:.4f}, risk {R[kb]:.2e} (b alone {cume[-1] / n:.2e});"
                         f" oracle per-eigvec keep/drop {orc:.2e}; M for target at best k: " +
                         (f"{tailv[kb] / (n * (TARGET - cume[kb] / n)):.0f}" if cume[kb] / n < TARGET else "unreachable (missed part alone > target)"))
    # alignment: where does the error sit in the variance spectrum?
    for bn, e in E.items():
        ec = (UA.T @ e) ** 2
        lines.append(f"   alignment b={bn}: error share in pilot eigvec bands [0,16) [16,64) [64,256) [256,512) [512,1024): " +
                     " ".join(f"{ec[a:b].sum() / ec.sum():.3f}" for a, b in ((0, 16), (16, 64), (64, 256), (256, 512), (512, n))) +
                     "; variance share: " + " ".join(f"{lamB_on_A[a:b].sum() / lamB_on_A.sum():.3f}" for a, b in ((0, 16), (16, 64), (64, 256), (256, 512), (512, n))))
    # coherent and weight-derived families
    T = tangents(W64, st)
    fam = {}
    a_L = W64[-1] @ (cl[-2] / np.linalg.norm(cl[-2]))
    for J in (2, 4, 8, 12): fam[f"coherent J={J}"] = hermite_feats(a_L, J)
    prop = [hermite_feats(a_L, 6)]
    for l in range(1, L - 1):
        a_l = W64[l] @ (cl[l - 1] / np.linalg.norm(cl[l - 1])); prop.append(T[l] @ hermite_feats(a_l, 6))
    fam["propagated coherent 15x6"] = np.concatenate(prop, 1)
    sdL = st[-1]["sd"]; aL = st[-1]["alpha"]; phL = np.exp(-0.5 * aL ** 2) / np.sqrt(2 * np.pi)
    He = [np.ones_like(aL), -aL]
    for j in range(2, 9): He.append(-aL * He[-1] - (j - 1) * He[-2])
    fam["readout shape 9"] = np.stack([sdL * phL * h for h in He], 1)
    fam["readout shape 9 + coherent 8"] = np.concatenate([fam["readout shape 9"], hermite_feats(a_L, 8)], 1)
    fam["CC1-closure direction"] = (cc[-1] - cl[-1])[:, None]
    for l in (0, 1, 3, 7, 11, 13):
        U, s, _ = np.linalg.svd(T[l]); e2 = s ** 2 / np.sum(s ** 2)
        for k in (8, 32, 128):
            fam[f"tangent T(L<-{l + 1}) top {k} (energy {e2[:k].sum():.2f})"] = U[:, :k]
    # local-defect-weighted Gram of all tangents: G = sum_l |eps_l|^2/n T_l T_l^T, eps_l from the truth (oracle sizes)
    G = np.zeros((n, n))
    for l in range(1, L):
        eps = (truth[l] - cl[l]) - (st[l]["Pa"] * (W64[l] @ (truth[l - 1] - cl[l - 1])))
        G += np.mean(eps ** 2) * (T[l] @ T[l].T)
    lg, Ug = np.linalg.eigh(G); Ug = Ug[:, ::-1]
    for k in (8, 32, 128, 256): fam[f"defect-weighted tangent Gram top {k}"] = Ug[:, :k]
    lines.append("(1b) captured fraction |P e|^2/|e|^2 per family (dim) and held-out variance tr(V^T S V)/n of the family:")
    for fn, X in fam.items():
        V = orth(X); vv = np.trace(V.T @ SB @ V) / n
        lines.append(f"   {fn:52s} r={V.shape[1]:4d}: " + "  ".join(f"{bn} {np.sum((V.T @ e) ** 2) / np.sum(e ** 2):.4f}" for bn, e in E.items()) +
                     f"   | var {vv:.2e} -> var/(M=4500) {vv / M_FLOOR:.1e}")
    for bn, e in E.items():
        v = e / np.linalg.norm(e); lines.append(f"   oracle rank-1 V = e/|e| (b={bn}): risk {v @ SB @ v / (n * M_FLOOR):.2e} (Rayleigh quotient {v @ SB @ v:.4f})")
    # ---------------- (2) projected residual covariance with exact-mean controls ----------------
    cE = math.exp(math.lgamma(d / 2) - math.lgamma((d + 1) / 2)) / (2 * math.sqrt(math.pi))   # E relu(theta_1)
    m1 = np.linalg.norm(W64[0], axis=1) * cE
    A1 = W64[0].T @ W64[0]; A1t = A1 - np.trace(A1) / d * np.eye(d)
    def feats(TH):
        return np.concatenate([TH.astype(np.float64), np.maximum(TH @ W32[0].T, 0).astype(np.float64) - m1,
                               np.sum((TH @ A1t.astype(np.float32)) * TH, 1, dtype=np.float64)[:, None]], 1)
    XA, XB = feats(THA), feats(THB); cols = {("lin",): d, ("lin", "h1"): d + n, ("lin", "h1", "gram"): d + n + 1}
    ctl = {"none": None}
    for which, p in cols.items():
        beta = np.linalg.solve(XA[:, :p].T @ XA[:, :p] + 1e-8 * np.eye(p), XA[:, :p].T @ (YA - mA).astype(np.float64))
        ctl["+".join(which)] = (beta, (YB - mB) - XB[:, :p] @ beta, p)
    lines.append("(2) held-out projected residual covariance tr(V^T S_z V)/n, and risk at M=4500, per control (fit on pilot):")
    SR = {}
    for cn, cv in ctl.items():
        SR[cn] = SB if cv is None else cov(cv[1], mean=cv[1].mean(0))[0]
    sel = [("all of R^n", np.eye(n)), ("coherent J=8", orth(fam["coherent J=8"])), ("propagated coherent 15x6", orth(fam["propagated coherent 15x6"])),
           ("defect-weighted tangent Gram top 256", fam["defect-weighted tangent Gram top 256"]),
           ("pilot quiet: bottom 512", UA[:, 512:]), ("pilot quiet: bottom 256", UA[:, 768:])]
    for nm, V in sel:
        row = []
        for cn, S_ in SR.items():
            v = np.trace(V.T @ S_ @ V) / n
            miss = {bn: (np.sum(e ** 2) - np.sum((V.T @ e) ** 2)) / n for bn, e in E.items()}
            row.append(f"{cn}: {v:.2e} (risk CC1 {miss['CC1'] + v / M_FLOOR:.1e})")
        lines.append(f"   {nm:40s} r={V.shape[1]:4d}  " + " | ".join(row))
    # oracle lower bounds (e known): best orthogonal projection of any rank (rank <= 1 since e e^T - S/M has one
    # positive eigenvalue at most) and best linear filter A (risk |e|^2 / (1 + M e^T S^-1 e)); raw and controlled
    lines.append("(1c) oracle lower bounds at M=4500 (e known exactly; no realizable estimator does better):")
    for cn in ("none", "lin+h1+gram"):
        S_ = SR[cn]; ls, Us = np.linalg.eigh(S_); keep = ls > 1e-10 * ls[-1]       # null directions (dead neurons) are exact for MC
        for bn, e in E.items():
            ve = e / np.linalg.norm(e); lmax = np.linalg.eigvalsh(np.outer(e, e) - S_ / M_FLOOR)[-1]
            ce = Us.T @ e; maha = np.sum(ce[keep] ** 2 / ls[keep]); nul = np.sum(ce[~keep] ** 2) / np.sum(ce ** 2)
            Rp = (e @ e - max(lmax, 0.0)) / n; Rl = np.sum(ce[keep] ** 2) / (1 + M_FLOOR * maha) / n; rq = ve @ S_ @ ve
            per = 32.8 * n * n / (2048 * n ** 3)                      # C/B per sample
            lines.append(f"   control {cn:12s} b={bn:8s}: best projection {Rp:.2e}, best linear filter {Rl:.2e} (M e^T S^+ e = {M_FLOOR * maha:.3f}; "
                         f"{(~keep).sum()} null directions hold {nul:.1e} of |e|^2), rank-1 V=e/|e| {rq / (n * M_FLOOR):.2e}"
                         f" (Rayleigh {rq:.3f} = {rq / (np.trace(S_) / n):.1f}x mean eigenvalue); adjusted score of the rank-1 oracle as M -> inf: {rq / n * per:.2e} (v56: 3.1e-09)")
    # quiet-subspace risk with the best control, re-optimised k (basis = pilot eigvecs of the controlled residual)
    bc = "lin+h1+gram"; RA = (YA - mA) - XA[:, :ctl[bc][2]] @ ctl[bc][0]; SRA = cov(RA, mean=RA.mean(0))[0]
    lr, Ur = np.linalg.eigh(SRA); Ur = Ur[:, ::-1]; lrB = np.sum(Ur * (SR[bc] @ Ur), 0)
    for bn, e in E.items():
        ec = (Ur.T @ e) ** 2; cume = np.concatenate([[0], np.cumsum(ec)]); tailv = np.concatenate([np.cumsum(lrB[::-1])[::-1], [0]])
        R = cume / n + tailv / (n * M_FLOOR); kb = int(np.argmin(R))
        lines.append(f"   controlled ({bc}) quiet subspace, b={bn}: best k={kb}, captured {1 - cume[kb] / cume[-1]:.4f}, risk {R[kb]:.2e} "
                     f"(b alone {cume[-1] / n:.2e}); variance reduction by the control on all of R^n: {np.trace(SB) / np.trace(SR[bc]):.2f}x")
    # ---------------- (3) metric codes: held-out int |z|^2 Psi_H and exact J(q) ----------------
    lines.append("(3) metric codes on held-out samples (z = V^T(Y - b) - c with V = R^n, b = CC1; w = |z|^2). A gain needs int w Psi_H / int w"
                 " of order d^2/16 (relative J reduction (4 G / (d E w))^2 / 8 at the best t); the ceiling bounds every density:")
    for cn in ("none", bc):
        Zc_A = (YA - bs["CC1"]) if cn == "none" else (RA + (mA - bs["CC1"]))
        Zc_B = (YB - bs["CC1"]) if cn == "none" else (ctl[cn][1] + (mB - bs["CC1"]))
        wA = np.einsum("ij,ij->i", Zc_A, Zc_A); wB = np.einsum("ij,ij->i", Zc_B, Zc_B)
        M0 = (THA.T * wA.astype(np.float32)) @ THA / len(wA) - (wA.mean() / d) * np.eye(d, dtype=np.float32)
        M0 = M0.astype(np.float64); M0 = 0.5 * (M0 + M0.T)
        # input-Jacobian Gram of the output, closure tangent from the input: A_in = T_(L<-0) diag(1/2) W_0 (d x d)
        Ain = T[-1]
        codes = {"W1^T W1": A1, "input-Jacobian Gram": Ain.T @ Ain, "pilot M0": M0}
        ve = E["CC1"] / np.linalg.norm(E["CC1"]); z1 = Zc_B @ ve; w1 = z1 * z1
        out = [f"   control {cn}: IS ceiling (any density q, q* ~ |z|): J(q*)/J(1) = (E|z|)^2/E|z|^2 = {np.mean(np.sqrt(wB)) ** 2 / wB.mean():.4f} for V = R^n,"
               f" {np.mean(np.abs(z1)) ** 2 / w1.mean():.4f} for the rank-1 oracle V = e_CC1/|e_CC1|",
               f"   control {cn}: E w = {wB.mean():.3e}; single ellipsoid along pilot M0: held-out first-order slope (d/2) tr(H M0_B)/E w = "]
        M0B = (THB.T * wB.astype(np.float32)) @ THB / len(wB) - (wB.mean() / d) * np.eye(d, dtype=np.float32)
        Hs = M0 / np.linalg.norm(M0); slope = (d / 2) * np.sum(Hs * M0B) / wB.mean()
        out[1] += f"{slope:+.3e} per unit t (|H|_F = 1)"
        hv, hU = np.linalg.eigh(Hs); pB = (THB.astype(np.float64) @ hU) ** 2
        jt = []
        for t in (-1.0, -0.5, -0.25, 0.25, 0.5, 1.0):
            qq = np.exp(-(d / 2) * np.log(pB @ np.exp(-t * hv)))          # det(e^{tH}) = e^{t tr H}, tr Hs ~ 0 (traceless M0)
            qq *= np.exp(-0.5 * t * hv.sum()); jt.append(f"t={t:+.1f}:{np.mean(wB / qq) / wB.mean():.4f}")
        out.append("      single-ellipsoid exact held-out J(q_t)/J(1): " + " ".join(jt))
        for hn, Mx in codes.items():
            ev, U = np.linalg.eigh(Mx); P = U[:, ev >= np.median(ev)][:, : d // 2]
            pa = np.sum((THB.astype(np.float64) @ P) ** 2, 1); s = 2 * pa - 1
            psi = (d * (d + 2) / 8) * s * s - d / 4
            g = wB * psi; G2, se = g.mean(), g.std() / math.sqrt(len(g))
            jt = []
            for tt in (0.01, 0.02, 0.04, 0.08):
                qm = 0.5 * (np.exp(-(d / 2) * np.log(np.exp(-tt) * pa + np.exp(tt) * (1 - pa))) + np.exp(-(d / 2) * np.log(np.exp(tt) * pa + np.exp(-tt) * (1 - pa))))
                jt.append(f"t={tt}:{np.mean(wB / qm) / wB.mean():.4f}")
            g1 = w1 * psi
            out.append(f"      H from {hn:22s}: int w Psi_H / int w = {G2 / wB.mean():+.3e} +- {se / wB.mean():.1e}; exact J(q_mix)/J(1) " + " ".join(jt) +
                       f" | rank-1 oracle z: {g1.mean() / w1.mean():+.2e} +- {g1.std() / math.sqrt(len(g1)) / w1.mean():.1e}")
        lines += out
    np.savez(f"{OUT}/s24_terminal_net{net}.npz", lamA=lamA, lamB_on_A=lamB_on_A, lamK=lamK, lamB_on_K=lamB_on_K,
             **{f"e_{k}": v for k, v in E.items()}, **{f"ecA_{k}": (UA.T @ v) ** 2 for k, v in E.items()},
             **{f"ecK_{k}": (UK.T @ v) ** 2 for k, v in E.items()}, mA=mA, mB=mB, trSR=np.array([np.trace(SR[c]) for c in SR]))
    txt = "\n".join(lines); print(txt, flush=True)
    with open(f"{OUT}/s24_terminal_net{net}.txt", "w") as f: f.write(txt + "\n")


if __name__ == "__main__":
    D, OUT = sys.argv[1], sys.argv[2]; V56 = os.environ.get("V56", "")
    for net in map(int, sys.argv[3:]): measure(D, OUT, net, V56)
