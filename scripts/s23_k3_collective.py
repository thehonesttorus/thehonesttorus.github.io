"""Is the per-neuron skewness of the final preactivations collective? (python scripts/s23_k3_collective.py KIK1 S21C DATA NET ...)

Exact final-row split Z_k = a_k c + zeta_k, a_k = w_k . mu_hat, c = mu_hat . h, zeta_k = w_perp_k . h_perp (E zeta_k = 0).
Averaging over the i.i.d. Gaussian perpendicular row parts (independent of the network below) gives the coherent cumulants
  E[k3(Z_k) | a_k] = a_k^3 k3(c) + 3 a_k s^2 Cov(c, r^2),          r^2 = |h|^2 - c^2,
  E[k4(Z_k) | a_k] = a_k^4 k4(c) + 6 a_k^2 s^2 [E(dc^2 r^2) - Var(c) E r^2] + 3 s^4 [Var(r^2) - 2|C_perp|_F^2] (the last from Wick),
so the dominant skewness is fixed by two collective scalars of the penultimate law. Compared with the Monte Carlo
per-neuron cumulants (kik1 pools; the collective moments from the s21c pools, all at the penultimate layer)."""
import sys, os, glob, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from scripts.kik_merge import pool_net
from scripts.s21_collective import MIX
rms = lambda x: np.sqrt(np.mean(np.square(x)))
K1, SC, D = sys.argv[1], sys.argv[2], sys.argv[3]
for net in map(int, sys.argv[4:]):
    a = pool_net(K1, net); N = a["N"]; P = a["pw_top"] / N
    t, q, m3, m4 = P[:, 0], P[:, 1], P[:, 2], P[:, 3]; var = q - t * t; sd = np.sqrt(var)
    k3 = m3 - 3 * t * q + 2 * t ** 3; k4 = m4 - 4 * t * m3 + 6 * t * t * q - 3 * t ** 4 - 3 * var ** 2
    g3, g4 = k3 / sd ** 3, k4 / var ** 2
    fs = sorted(glob.glob(f"{SC}/s21c_{net}_*.npz")); NN = 0; mix = 0
    for f in fs: z = np.load(f); NN += int(z["N"]); mix = mix + z["mix"]
    M = {ij: mix[14, kk] / NN for kk, ij in enumerate(MIX)}           # penultimate post-activation (0-based layer 14)
    Ec, Ec2, Ec3, Ec4 = M[(1, 0)], M[(2, 0)], M[(3, 0)], M[(4, 0)]; Er, Er2, Ecr, Ec2r = M[(0, 1)], M[(0, 2)], M[(1, 1)], M[(2, 1)]
    Vc = Ec2 - Ec ** 2; k3c = Ec3 - 3 * Ec * Ec2 + 2 * Ec ** 3; k4c = Ec4 - 4 * Ec * Ec3 + 6 * Ec * Ec * Ec2 - 3 * Ec ** 4 - 3 * Vc ** 2
    cov_cr = Ecr - Ec * Er; Edc2r = Ec2r - 2 * Ec * Ecr + Ec * Ec * Er; Vr = Er2 - Er ** 2
    truth = np.load(f"{D}/truth_off{net}.npz")["m"].astype(np.float64); W = np.load(f"{D}/W_off{net}.npy")[-1].astype(np.float64)
    mh = truth[-2] / np.linalg.norm(truth[-2]); ak = W @ mh; n = len(ak); s2 = 2.0 / n
    k3_coh = ak ** 3 * k3c + 3 * ak * s2 * cov_cr
    g3_pred = k3_coh / sd ** 3
    g3_pred_cov = 3 * ak * s2 * cov_cr / sd ** 3
    noise3 = np.sqrt(6.0 / N)
    def share(y, x):
        r = y - x; return 1 - np.var(r) / np.var(y), np.polyfit(x, y, 1)[0]
    sh, sl = share(g3, g3_pred)
    print(f"net {net}: penultimate collective: Var c / Ec^2 {Vc / Ec ** 2:.4f}, skew c {k3c / Vc ** 1.5:+.3f}, corr(c, r^2) {cov_cr / np.sqrt(Vc * Vr):+.3f}")
    print(f"   skewness gamma3_k: rms {rms(g3):.4f} (MC noise {noise3:.4f}); predicted coherent rms {rms(g3_pred):.4f} (Cov(c,r^2) term alone {rms(g3_pred_cov):.4f})")
    print(f"   variance of gamma3 explained by the two-scalar collective prediction: {sh:.3f}; regression slope {sl:.3f}; residual rms {rms(g3 - g3_pred):.4f}")
    # kurtosis: coherent prediction needs |C_perp|_F^2; use the measured Var(r^2) split by the Wick formula at the ensemble level
    k4_mean_pred_scale = 3 * s2 ** 2 * Vr                       # upper bound without the Gaussian-part subtraction
    print(f"   excess kurtosis: mean {g4.mean():+.4f}; naive scale-mixture 3 s^4 Var(r^2)/sigma^4 mean {np.mean(k4_mean_pred_scale / var ** 2):+.4f} (the Wick term -6 s^4 |C_perp|_F^2 must be subtracted)")
