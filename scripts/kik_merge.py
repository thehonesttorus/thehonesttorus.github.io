"""Kikuchi certificate pipeline on the real networks, from the pooled single-pass contractions of kik_contract.py.

  python scripts/kik_merge.py RESULTS DATA NET [NET ...]      prints the report, writes RESULTS/kikm_{NET}.npz

For every final coordinate j (the same quantities at every layer for the spectral part):
  first interval [t_+, (t+s)/2] (Thm 4.1), restart effect E_nu g_c (Thm 5.1), quartic gap d4/(4s^3) (eq. 11),
  even/odd cyclic hierarchy L_k <= m <= U_k with residual R_k, k = 1..4 at the last layer, 1..2 below (Thms 7.1-8.1),
  nested conditional cells on K penultimate gates (Thm 10.1, eqs. 21-22) and the best single-gate refinement,
  pair overlap m = E P - kappa E_nu 1/max(P, Q) (Prop. 9.1),
  and the same first interval from the Gaussian-closure contractions (the analytic first pass).
Widths are reported as the midpoint RMS certificate ||u - l|| / (2 sqrt n) (eq. 21); the target RMS is 6e-5.
All Monte Carlo contractions are exact moments of the pooled empirical law, so each certificate is exact for that
law; "err" columns compare midpoints with the true means of the truth file (whose own MC error is reported too).
"""
import sys, os, glob, numpy as np
from multiprocessing import Pool
from scipy.special import ndtr
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from whest import kikuchi as kk


def pool_net(R, net):
    fs = sorted(glob.glob(f"{R}/kc_{net}_*.npz")); z0 = np.load(fs[0]); acc = {k: z0[k].astype(np.float64) for k in z0.files if k not in ("order",)}
    order = z0["order"]
    for f in fs[1:]:
        z = np.load(f); assert np.array_equal(z["order"], order), "pilot orders differ"
        for k in acc: acc[k] = acc[k] + z[k]
    acc["order"] = order; acc["ntask"] = len(fs); return acc


def gauss_pos(t, q):
    sd = np.sqrt(np.maximum(q - t * t, 1e-300)); a = t / sd
    return t * ndtr(a) + sd * np.exp(-0.5 * a * a) / np.sqrt(2 * np.pi)


def _spec(args):
    t, M, kmax = args
    return [kk.spectral(t, M, k, dps=40) for k in range(1, kmax + 1)]


def rmsw(lo, hi): return kk.rms_bound(lo, hi)
def rms(x): return np.sqrt(np.mean(x ** 2))


def analyse(R, D, net, procs=4):
    a = pool_net(R, net); N = a["N"]; truth = np.load(f"{D}/truth_off{net}.npz")["m"].astype(np.float64); L, n = truth.shape
    mc = a["pos"] / N; out = {}
    pw = [np.concatenate([np.ones((n, 1)), a["pw_low"][l] / N], 1) for l in range(L - 1)] + [np.concatenate([np.ones((n, 1)), a["pw_top"] / N], 1)]
    t, q = pw[-1][:, 1], pw[-1][:, 2]; m = mc[-1]; s = np.sqrt(q); tr = truth[-1]
    print(f"\n=== net {net}: {N} samples in {a['ntask']} tasks; final layer: mean t {t.mean():+.3f}, rms t {rms(t):.3f}, mean s {s.mean():.3f}")
    print(f"  Monte Carlo m vs truth: rms {rms(m - tr):.2e} (sampling sd ~ {np.sqrt(np.mean(q) / N):.1e} upper bound)")
    # Thm 4.1
    lo, hi, Dl = kk.moment_interval(t, q)
    print(f"  [4.1] first interval: RMS cert {rmsw(lo, hi):.3e}; mean half-gap Delta {Dl.mean():.3f}; midpoint err {rms((lo + hi) / 2 - tr):.3e}")
    # Thm 5.1 restart effect
    eff = (m - lo) / Dl
    print(f"  [5.1] restart: E_nu g_c = (m - t_+)/Delta in [{eff.min():.3f}, {eff.max():.3f}], mean {eff.mean():.3f}, sd {eff.std():.3f}"
          f" (a certified [l,u] on it gives RMS cert Delta*(u-l)/2: needs u-l <= {1.2e-4 / rms(Dl):.1e})")
    # eq. 11
    r4 = pw[-1][:, 4]; Qd = kk.quartic(q, r4); gap = hi - m
    print(f"  [6]   quartic: d4/q^2 mean {np.mean(Qd['d4'] / q ** 2):.3f} (Gaussian with these t,q: {np.mean((3 * (q - t * t) ** 2 + 6 * (q - t * t) * t * t + t ** 4 - q * q) / q ** 2):.3f});"
          f" gap {gap.mean():.4f} <= bound {Qd['gap_bound'].mean():.4f}; E_nu4 effect mean {np.mean(gap * 4 * s ** 3 / Qd['d4']):.3f}")
    # Gaussian shape at the exact first two moments (diagnostic) and the analytic first pass
    gp = gauss_pos(t, q)
    print(f"  [ref] Gaussian shape at exact (t, q): err vs truth {rms(gp - tr):.3e}; vs the same-sample m {rms(gp - m):.3e} (the shape defect itself)")
    # how much shape information the next cumulants carry (same-sample, so sampling noise largely cancels)
    sd = np.sqrt(q - t * t); al = t / sd; ph = np.exp(-al * al / 2) / np.sqrt(2 * np.pi)
    c3 = (pw[-1][:, 3] - 3 * t * q + 2 * t ** 3) / sd ** 3
    mu4 = pw[-1][:, 4] - 4 * t * pw[-1][:, 3] + 6 * t * t * q - 3 * t ** 4; c4 = mu4 / sd ** 4 - 3
    d0 = gauss_pos(t, q) / sd; d3 = -al * ph; d4 = (al * al - 1) * ph; d6 = (al ** 4 - 6 * al * al + 3) * ph
    e3 = sd * (d0 + c3 / 6 * d3); e34 = sd * (d0 + c3 / 6 * d3 + c4 / 24 * d4 + c3 ** 2 / 72 * d6)
    eff_g = (gp - lo) / Dl; eff = (m - lo) / Dl
    print(f"        standardised kappa3 rms {rms(c3):.3f}, kappa4 mean {c4.mean():+.3f}; Hermite shape with kappa3: err {rms(e3 - m):.3e}, kappa3+kappa4: {rms(e34 - m):.3e};"
          f" restart effect: Gaussian value err {rms(eff_g - eff):.2e} (needs <= {1.2e-4 / rms(Dl):.1e} for the target)")
    out.update(c3=c3, c4=c4, gp=gp)
    gfile = f"{R}/kcg_{net}.npz"
    if os.path.exists(gfile):
        g = np.load(gfile); tg, qg = g["t"][-1], g["q"][-1]; lg, hg, _ = kk.moment_interval(tg, qg)
        cover = np.mean((tr >= lg - 1e-12) & (tr <= hg + 1e-12))
        print(f"  [A2]  closure contractions: t err {rms(tg - t):.2e}, q err {rms(qg - q):.2e}; interval RMS cert {rmsw(lg, hg):.3e}, covers truth {cover:.3f};"
              f" closure mean err {rms(g['m'][-1] - tr):.3e}; Gaussian at closure (t,q) err {rms(gauss_pos(tg, qg) - tr):.3e}")
        out.update(tg=tg, qg=qg)
    # spectral hierarchy, every layer
    jobs = [(pw[l][j, 1], pw[l][j], 4 if l == L - 1 else 2) for l in range(L) for j in range(n)]
    with Pool(procs) as P:
        res = P.map(_spec, jobs, chunksize=64)
    print("  [7-8] cyclic hierarchy (RMS certs; maintained interval [max(t_+, L_k, U_k - R_k), U_k] intersected over k):")
    for l in list(range(0, L - 1, 3)) + [L - 1]:
        rs = res[l * n:(l + 1) * n]; kmax = len(rs[0]); line = []; LO = np.maximum(pw[l][:, 1], 0); HI = (pw[l][:, 1] + np.sqrt(pw[l][:, 2])) / 2
        for k in range(kmax):
            U = np.array([r[k]["U"] for r in rs]); Lk = np.array([r[k]["L"] for r in rs]); Rk = np.array([r[k]["R"] for r in rs])
            LO = np.maximum(LO, np.maximum(Lk, U - Rk)); HI = np.minimum(HI, U)
            line.append(f"k={k + 1}: [L,U] {rmsw(Lk, U):.2e} kept {rmsw(LO, HI):.2e} mid err {rms((LO + HI) / 2 - truth[l]):.2e} R/(U-L) {np.median(Rk / (U - Lk)):.2f}")
            if l == L - 1: out[f"spec_lo{k + 1}"], out[f"spec_hi{k + 1}"] = LO.copy(), HI.copy()
        print(f"    layer {l:2d}: " + " | ".join(line))
    # Thm 10.1 nested cells on the K pilot-ordered penultimate gates
    K = a["order"].shape[1]; cnt = a["cnt"].reshape(n, 1 << K) / N; cz = a["cz"].reshape(n, 1 << K) / N
    cz2 = a["cz2"].reshape(n, 1 << K) / N; czp = a["czp"].reshape(n, 1 << K) / N
    print(f"  [10]  nested cells on {K} penultimate gates per final neuron (pilot greedy order):")
    lines = []
    for k in range(K + 1):
        p_ = cnt.reshape(n, -1, 1 << k).sum(1); t_ = cz.reshape(n, -1, 1 << k).sum(1); q_ = cz2.reshape(n, -1, 1 << k).sum(1); zp = czp.reshape(n, -1, 1 << k).sum(1)
        lo_k, hi_k, _ = kk.cell_interval(t, p_, t_, q_)
        resid = np.abs((m - lo_k) - np.sum(np.minimum(zp, zp - t_), 1)).max()
        lines.append(f"K={k}: cert {rmsw(lo_k, hi_k):.3e} (x{rmsw(lo_k, hi_k) / rmsw(lo, hi):.3f}) mid err {rms((lo_k + hi_k) / 2 - tr):.2e}")
        out[f"cell_lo{k}"], out[f"cell_hi{k}"] = lo_k, hi_k
    print("    " + "\n    ".join(lines)); print(f"    exact lower residual identity, max error {resid:.1e}")
    M0 = a["M0"] / N; M1 = a["M1"] / N; M2 = a["M2"] / N
    gw = lambda p_, tt, qq: 0.5 * (np.sqrt(np.maximum(p_ * qq, 0)) - np.abs(tt))
    dec = gw(1.0, t[None], q[None]) - gw(M0[:, None], M1, M2) - gw(1 - M0[:, None], t[None] - M1, q[None] - M2)
    best = dec.max(0); w1 = gw(1.0, t, q)
    print(f"    best single penultimate gate (full-sample eq. 22 scores): width decrease {np.median(best / w1):.4f} of the first width (median), max {np.max(best / w1):.4f}")
    # Prop 9.1 pair overlap
    EP, EQ, EP2, EQ2, EPQ, Emin, Eratio = a["pair"] / N
    corr = EPQ - EP * EQ
    print(f"  [9]   pair overlap: E P {EP.mean():.2f}, E Q {EQ.mean():.2f}, kappa {EPQ.mean():.1f}, corr(P,Q) {np.mean(corr / np.sqrt((EP2 - EP ** 2) * (EQ2 - EQ ** 2))):.3f};"
          f" identity err {np.abs(m - (EP - Emin)).max():.1e}; E_nu 1/max {np.mean(Emin / EPQ):.4f} vs 1/E max {np.mean(1 / (EP + EQ - Emin)):.4f};"
          f" relative spread of max under nu (from E_nu 1/max^2 vs (E_nu 1/max)^2): "
          f"{np.mean(np.sqrt(np.maximum(Eratio / EPQ - (Emin / EPQ) ** 2, 0)) / (Emin / EPQ)):.4f}")
    np.savez(f"{R}/kikm_{net}.npz", t=t, q=q, m=m, truth=tr, **out)


if __name__ == "__main__":
    R, D = sys.argv[1], sys.argv[2]
    for net in map(int, sys.argv[3:]): analyse(R, D, net)
