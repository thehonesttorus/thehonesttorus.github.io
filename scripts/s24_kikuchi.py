"""The Kikuchi tower of a real network: levels, level transitions, arrow weights and gate-law temperature, per layer.

  python scripts/s24_kikuchi.py DATA OUT NET [NSAMP]          writes OUT/s24_kikuchi_net{NET}.txt and .npz

x ~ N(0, I) and its Mehler rotation partner x_t = x cos t + y sin t (t = T_ROT) are pushed through the network.
Per layer l (gates g_l = 1{z_l > 0}, level R_l = |S_l| = sum_i g_(l,i)):
  levels       mean and sd of R_l / n; corr(R_l, R_(l+1)); the share of Var R_(l+1) explained linearly by the signs
               g_l, by the values a_l, and by both (the next level depends on values, not only on signs);
  arrows       flip counts of every gate between x and x_t: flip rate r_(l,i) = flips / (N t) (both directions; Rice's
               formula: r = E[|zdot| delta(z)]), up/down balance, and E zdot^2 (= E |grad z|^2) against E z^2;
               the Gaussian Rice prediction (1/pi) exp(-alpha^2/2) sqrt(E zdot^2)/sd with the closure's alpha, sd;
  temperature  the gate correlation matrix R_g: its top eigenvalue (eta = lambda_max - 1 is the ALO influence radius
               of the unpinned law) and participation ratio, and eta after regressing the gates on the top-k principal
               components of the pre-activations (k collective coordinates conditioned out); the MP noise edge for
               N samples is (1 + sqrt(n/N))^2 - 1;
  depth        the closure's graded gap g_l = 2 E Phi(alpha)^2 ... (Stage 21) and the depth gain prod_(k > l) g_k."""
import sys, os, time, math, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from whest import gate_channel as gc

BS, T_ROT, KS = 4096, 0.01, (0, 1, 2, 4, 8, 16, 32, 64, 128)
from scipy.special import ndtr
ndtr_abs = lambda a: ndtr(np.abs(a))


def run(D, OUT, net, N):
    W = [np.ascontiguousarray(w, dtype=np.float32) for w in np.load(f"{D}/W_off{net}.npy")]; L, n = len(W), W[0].shape[0]
    d = W[0].shape[1]
    Sg = np.zeros((L, n)); Sgg = np.zeros((L, n, n)); Sz = np.zeros((L, n)); Szz = np.zeros((L, n, n)); Sgz = np.zeros((L, n, n))
    Sa = np.zeros((L, n)); Saa = np.zeros((L, n, n)); Sga = np.zeros((L, n, n))
    SR = np.zeros(L); SRR = np.zeros(L); SRnext = np.zeros(L); SgRn = np.zeros((L, n)); SaRn = np.zeros((L, n))
    up = np.zeros((L, n)); dn = np.zeros((L, n)); szd2 = np.zeros((L, n)); sz2 = np.zeros((L, n)); szdf = np.zeros((L, n))
    DEL = (1, 2, 4); Sx = {dl: np.zeros((L, n, n)) for dl in DEL}
    SRt = np.zeros(L); SE = np.zeros(L); SEE = np.zeros(L); SdR2 = np.zeros(L)
    lvb = np.linspace(-2.5, 2.5, 6); lup = np.zeros((L, 7)); ldn = np.zeros((L, 7)); locc = np.zeros((L, 7)); lact = np.zeros((L, 7))
    st = gc.closure_states([w.astype(np.float64) for w in W], keep_C=False)
    ref = [s_["alpha"] > 0 for s_ in st]                                     # reference (majority) pattern of the closure
    mu_lvl = np.array([np.sum(s_["Pa"]) for s_ in st]); sd_lvl = np.full(L, 1.0)
    t0 = time.time()
    for b in range(N // BS):
        r = np.random.default_rng([31_000 + net, b])
        X = r.standard_normal((d, BS), dtype=np.float32); Y = r.standard_normal((d, BS), dtype=np.float32)
        Xt = X * np.float32(math.cos(T_ROT)) + Y * np.float32(math.sin(T_ROT))
        h, ht, Rs, gs, As = X, Xt, [], [], []
        for l, Wl in enumerate(W):
            z = Wl @ h; zt = Wl @ ht; g = (z > 0); gt = (zt > 0); h = np.maximum(z, 0); ht = np.maximum(zt, 0)
            gf = g.astype(np.float32); R = gf.sum(0, dtype=np.float64)
            Sg[l] += gf.sum(1); Sgg[l] += gf @ gf.T; Sz[l] += z.sum(1); Szz[l] += z @ z.T; Sgz[l] += gf @ z.T
            Sa[l] += h.sum(1); Saa[l] += h @ h.T; Sga[l] += gf @ h.T
            SR[l] += R.sum(); SRR[l] += R @ R
            up[l] += np.sum(gt & ~g, 1); dn[l] += np.sum(g & ~gt, 1)
            zd = (zt - z) / T_ROT; szd2[l] += np.sum(zd.astype(np.float64) ** 2, 1); sz2[l] += np.sum(z.astype(np.float64) ** 2, 1)
            szdf[l] += np.sum(np.abs(zd) * (g != gt), 1)
            Rs.append(R); gs.append(gf); As.append(h)
            Rt = gt.sum(0, dtype=np.float64); SRt[l] += R @ Rt; SdR2[l] += np.sum((R - Rt) ** 2)
            ex = np.sum(g != ref[l][:, None], 0, dtype=np.float64); SE[l] += ex.sum(); SEE[l] += ex @ ex
            if b == 0: sd_lvl[l] = max(R.std(), 1.0)
            bi = np.digitize((R - mu_lvl[l]) / sd_lvl[l], lvb)                 # level bins in units of the level sd
            for k in range(7):
                m_ = bi == k
                if m_.any():
                    locc[l, k] += m_.sum(); lact[l, k] += R[m_].sum()
                    lup[l, k] += np.sum((gt & ~g)[:, m_]); ldn[l, k] += np.sum((g & ~gt)[:, m_])
        for l in range(L - 1):
            SRnext[l] += Rs[l] @ Rs[l + 1]; SgRn[l] += gs[l] @ Rs[l + 1]; SaRn[l] += As[l] @ Rs[l + 1]
        for dl in DEL:
            for l in range(L - dl): Sx[dl][l] += gs[l] @ gs[l + dl].T
    tmc = time.time() - t0
    lines = [f"\n=== net {net}: N = {N} Gaussian inputs and rotation partners (t = {T_ROT}), {tmc:.0f}s ===",
             " l | level R/n: mean  sd   | corr(R_l,R_l+1) | R^2 of R_l+1 on signs / values / both | flip rate: mean  sd/mean  up/down"
             " | E zdot^2/E z^2 | Rice-Gauss pred./measured | gate corr: lmax-1  PR | eta after k PCs out (k=" + ",".join(map(str, KS)) + ")"
             " | MP edge | gap g_l  prod_(k>l) g_k"]
    gaps = np.array([2 * np.mean(s["Pa"] ** 2) for s in st]); out = {}
    mp = (1 + math.sqrt(n / N)) ** 2 - 1
    for l in range(L):
        mg = Sg[l] / N; Cg = Sgg[l] / N - np.outer(mg, mg); mz = Sz[l] / N; Cz = Szz[l] / N - np.outer(mz, mz)
        Cgz = Sgz[l] / N - np.outer(mg, mz)
        mR = SR[l] / N; vR = SRR[l] / N - mR * mR
        if l < L - 1:
            mR1 = SR[l + 1] / N; vR1 = SRR[l + 1] / N - mR1 ** 2; cRR = (SRnext[l] / N - mR * mR1) / math.sqrt(vR * vR1)
            ma = Sa[l] / N; Ca = Saa[l] / N - np.outer(ma, ma); Cga = Sga[l] / N - np.outer(mg, ma)
            cg = SgRn[l] / N - mg * mR1; ca = SaRn[l] / N - ma * mR1
            r2 = lambda C, c: float(c @ np.linalg.lstsq(C + 1e-9 * np.trace(C) / n * np.eye(len(c)), c, rcond=None)[0] / vR1)
            Cb = np.block([[Cg, Cga], [Cga.T, Ca]]); cb = np.concatenate([cg, ca])
            lev = f"{cRR:+.3f} | {r2(Cg, cg):.3f} / {r2(Ca, ca):.3f} / {r2(Cb, cb):.3f}"
        else:
            lev = "   -   |   -   /   -   /   -  "
        rate = (up[l] + dn[l]) / (N * T_ROT); ezd2 = szd2[l] / N; ez2 = sz2[l] / N
        a, sd = st[l]["alpha"], st[l]["sd"]; pred = np.exp(-0.5 * a * a) / np.pi * np.sqrt(ezd2) / sd
        keep = rate > 0
        sg = np.sqrt(np.clip(np.diag(Cg), 1e-30, None)); live = np.diag(Cg) > 1e-6
        Rg = Cg[np.ix_(live, live)] / np.outer(sg[live], sg[live]); ev = np.linalg.eigvalsh(Rg)
        pr = float(np.sum(ev) ** 2 / np.sum(ev ** 2))
        lz, Uz = np.linalg.eigh(Cz); Uz = Uz[:, ::-1]
        etas = []
        for k in KS:
            if k == 0: Cr = Cg
            else:
                U = Uz[:, :k]; A = Cgz @ U; Cr = Cg - A @ np.linalg.solve(U.T @ Cz @ U, A.T)
            Cr = Cr[np.ix_(live, live)]; sr = np.sqrt(np.clip(np.diag(Cr), 1e-30, None))
            etas.append(float(np.linalg.eigvalsh(Cr / np.outer(sr, sr))[-1] - 1))
        lv = live & (rate > 0); Gm = rate[lv] / 2; Mg = Cg[np.ix_(lv, lv)] / np.sqrt(np.outer(Gm, Gm)); lin_gap = 1 / np.linalg.eigvalsh(Mg)[-1]
        lin_gap_k = []
        for k in (16, 128):
            U = Uz[:, :k]; A = Cgz @ U; Cr = (Cg - A @ np.linalg.solve(U.T @ Cz @ U, A.T))[np.ix_(lv, lv)]
            lin_gap_k.append(1 / np.linalg.eigvalsh(Cr / np.sqrt(np.outer(Gm, Gm)))[-1])
        out[f"lingap_{l}"] = np.array([lin_gap] + lin_gap_k)
        lines.append(f"{l + 1:2d} | {mR / n:.3f} {math.sqrt(vR) / n:.4f} | {lev} | {rate.mean():.4f} {rate.std() / rate.mean():.3f} "
                     f"{up[l].sum() / max(dn[l].sum(), 1):.3f} | {np.mean(ezd2) / np.mean(ez2):.3f} | "
                     f"{np.sum(pred) / np.sum(rate):.3f} | {ev[-1] - 1:7.2f} {pr:6.1f} | " + " ".join(f"{e:.2f}" for e in etas) +
                     f" | {mp:.2f} | {gaps[l]:.3f} {np.prod(gaps[l + 1:]):.3f} | wall-chain level-1 gap {lin_gap:.3f} (16/128 PCs out: {lin_gap_k[0]:.3f} {lin_gap_k[1]:.3f})")
        out[f"rate_{l}"] = rate; out[f"pred_{l}"] = pred; out[f"eta_{l}"] = np.array(etas); out[f"ev_{l}"] = ev
        out[f"ezd2_{l}"] = ezd2; out[f"ez2_{l}"] = ez2
    lines.append(" l | excitation |S - Sref|/n: mean sd (closure sum Phi(-|alpha|)/n) | level relaxation rate (1-corr(R(x),R(x_t)))/t"
                 " vs independent 2/pi | per-neuron rates by level bin (-2.5..2.5 sd): up per inactive / down per active"
                 " | cross-layer gate corr ||Corr(g_l, g_l+d)||_op d=1,2,4 (raw; MP edge sqrt(n/N)*2)")
    for l in range(L):
        mR = SR[l] / N; vR = SRR[l] / N - mR * mR; cR = 1 - SdR2[l] / N / (2 * vR)                 # 1 - corr = E(R - R_t)^2 / (2 Var R)
        mE = SE[l] / N; sE = math.sqrt(max(SEE[l] / N - mE * mE, 0)); pE = np.sum(1 - ndtr_abs(st[l]["alpha"])) / n
        rates = []
        for k in range(1, 6):
            if locc[l, k] > 50:
                ina = locc[l, k] * n - lact[l, k]; rates.append(f"{lup[l, k] / (ina * T_ROT):.3f}/{ldn[l, k] / (lact[l, k] * T_ROT):.3f}")
            else: rates.append("  -  ")
        mg = Sg[l] / N; sg = np.sqrt(np.clip(Sgg[l].diagonal() / N - mg * mg, 1e-30, None)); xs = []
        for dl in DEL:
            if l + dl < L:
                mg2 = Sg[l + dl] / N; sg2 = np.sqrt(np.clip(Sgg[l + dl].diagonal() / N - mg2 * mg2, 1e-30, None))
                Cx = (Sx[dl][l] / N - np.outer(mg, mg2)) / np.outer(sg, sg2); live = (sg > 1e-3)[:, None] & (sg2 > 1e-3)[None, :]
                xs.append(f"{np.linalg.norm(np.where(live, Cx, 0), 2):.2f}")
            else: xs.append("  - ")
        lines.append(f"{l + 1:2d} | {mE / n:.4f} {sE / n:.4f} ({pE:.4f}) | {(1 - cR) / T_ROT:6.2f} vs {2 / np.pi:.3f} | " + " ".join(rates) + " | " + " ".join(xs) +
                     f" (edge {2 * math.sqrt(n / N):.2f})")
    txt = "\n".join(lines); print(txt, flush=True)
    with open(f"{OUT}/s24_kikuchi_net{net}.txt", "w") as f: f.write(txt + "\n")
    np.savez(f"{OUT}/s24_kikuchi_net{net}.npz", gaps=gaps, **out)


if __name__ == "__main__":
    D, OUT, net = sys.argv[1], sys.argv[2], int(sys.argv[3]); N = int(sys.argv[4]) if len(sys.argv) > 4 else 1 << 17
    run(D, OUT, net, N)
