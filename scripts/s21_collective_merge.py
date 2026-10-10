"""Pool s21_collective.py tasks: python scripts/s21_collective_merge.py RESULTS DATA NET [NET ...]

Per selected layer l: the exact coherent function Ybar_l(a_j) = E_x G(a_j c, s^2 r^2) against the true means
(residual eta = incoherent part, its RMS, and the Wick-orthogonality statistics T_psi, which must be N(0, 1));
the size of the collective fluctuation (Ybar against G(a E c, s^2 E r^2)); the collective law of (c, r^2) per layer
(relative variances, correlation, standardised third cumulants); and the closure's coherent error: the closure means
against Ybar, and the closure's own collective scalars E c = |mu|, E r^2 = Tr S^u - mu_hat S^u mu_hat."""
import sys, os, glob, numpy as np
from scipy.special import ndtr
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from scripts.s21_collective import G, MIX, YB_LAYERS

rms = lambda x: np.sqrt(np.mean(np.square(x)))
R, D = sys.argv[1], sys.argv[2]
for net in map(int, sys.argv[3:]):
    fs = sorted(glob.glob(f"{R}/s21c_{net}_*.npz")); z0 = np.load(fs[0]); N = 0; mix = 0; yb = {l: 0 for l in YB_LAYERS}
    for f in fs:
        z = np.load(f); N += int(z["N"]); mix = mix + z["mix"]
        for l in YB_LAYERS: yb[l] = yb[l] + z[f"yb_{l}"]
    truth = np.load(f"{D}/truth_off{net}.npz")["m"].astype(np.float64); L, n = truth.shape; s2 = 2.0 / n
    M = {ij: mix[:, k] / N for k, ij in enumerate(MIX)}
    g = np.load(f"{R}/kcg_{net}.npz") if os.path.exists(f"{R}/kcg_{net}.npz") else None
    print(f"\n=== net {net}: {N} samples, {len(fs)} tasks")
    print(" layer | rms eta (incoherent) | T_psi z (1, He1, He2) | collective effect rms | closure err: total, coherent part")
    for l in YB_LAYERS:
        Y = yb[l] / N; a = z0[f"a_{l}"]; eta = truth[l] - Y; x = a / a.std()
        T = [np.mean(eta * p) / (np.std(eta * p) / np.sqrt(n)) for p in (np.ones(n), x, x * x - 1)]
        Ec, Er2 = M[(1, 0)][l - 1], M[(0, 1)][l - 1]; Y0 = G(a * Ec, s2 * Er2)
        line = f"  {l:2d}  | {rms(eta):.3e} | {T[0]:+.2f} {T[1]:+.2f} {T[2]:+.2f} | {rms(Y - Y0):.3e}"
        if g is not None:
            ec = g["m"][l] - truth[l]; coh = np.polyval(np.polyfit(x, ec, 8), x)   # E[closure - truth | a], by regression
            line += f" | {rms(ec):.2e}, {rms(coh):.2e}"
        print(line)
    print(" collective law per layer (post-activation l): E c = |mu| check, rel var c, rel var r^2, corr(c, r^2), skew c, skew r^2 | closure: E c, E r^2 rel err")
    for l in (0, 1, 2, 4, 8, 12, 14, 15):
        Ec, Ec2, Ec3 = M[(1, 0)][l], M[(2, 0)][l], M[(3, 0)][l]; Er, Er2_, Er3 = M[(0, 1)][l], M[(0, 2)][l], M[(0, 3)][l]
        vc = Ec2 - Ec ** 2; vr = Er2_ - Er ** 2; cov = M[(1, 1)][l] - Ec * Er
        sk = lambda m1, m2, m3: (m3 - 3 * m1 * m2 + 2 * m1 ** 3) / (m2 - m1 ** 2) ** 1.5
        line = (f"  {l:2d}: {Ec / np.linalg.norm(truth[l]) - 1:+.1e}  {vc / Ec ** 2:.4f}  {vr / Er ** 2:.4f}  {cov / np.sqrt(vc * vr):+.3f}  "
                f"{sk(Ec, Ec2, Ec3):+.3f}  {sk(Er, Er2_, Er3):+.3f}")
        print(line)
