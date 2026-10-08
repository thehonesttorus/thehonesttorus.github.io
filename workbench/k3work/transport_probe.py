# Variance of the Gaussian-transport estimator of a cumulant source's effect on the output, at real magnitudes.
#   python transport_probe.py NET MCPREFIX CUT [CUT ...]  [--n N]
# At the pre-activation z_b of layer b take the Gaussian reference N(mu_b, Sigma_b) (Monte Carlo truth, full covariance)
# and the true kappa3 and kappa4 diagonals as a source. Its first-order effect on the network's output means is
#   E[h F] with h = sum_a k3_a He3_a / 6 + k4_a He4_a / 24   (score estimator), and, by one Gaussian integration by parts,
#   E[U . grad F] with U_a = k3_a (v_a^2 - P_aa)/6 + k4_a (v_a^3 - 3 P_aa v_a)/24, v = P (z - mu), P = Sigma^-1
# (transport estimator: one tangent channel through the actual suffix, actual gates). Antithetic pairs (w, -w).
# Reports, per cut: the rms over output neurons of the signal (the transport mean), of each estimator's per-sample
# standard deviation, the samples needed for 10% relative error at the rms level, and the FLOPs that would buy.
import sys, math, time, numpy as np
net, pre = int(sys.argv[1]), sys.argv[2]
args = sys.argv[3:]; N = 4096
if "--n" in args:
    N = int(args[args.index("--n") + 1]); args = args[:args.index("--n")]
cuts = [int(c) for c in args]
Wcol = np.load(f"../official/W_off{net}.npy").astype(np.float64)
L, n, _ = Wcol.shape
F = np.load(f"{pre}_full.npz")
rng = np.random.default_rng(11)
for b in cuts:
    t0 = time.time()
    mu = F["mu"][b].astype(np.float64); S = F["cov"][b].astype(np.float64); S = 0.5 * (S + S.T)
    k3 = F["k3"][b].astype(np.float64); k4 = F["k4"][b].astype(np.float64)
    # the covariance is singular (dead units): restrict the reference to its support, z = mu + Q_r Lam_r^(1/2) w,
    # P = pseudo-inverse, and project the Stein field onto the range (the true kappa3 of z_b lies in it)
    lam, Q = np.linalg.eigh(S)
    keep = lam > 1e-6 * lam.max()
    Qr, lr = Q[:, keep], lam[keep]
    Lc = Qr * np.sqrt(lr)[None, :]                # z - mu = Lc w
    Linv = (Qr / np.sqrt(lr)[None, :]).T          # v = P (z - mu) = Linv^T w
    Pd = np.sum((Qr / np.sqrt(lr)[None, :]) ** 2, axis=1)
    Pi = Qr @ Qr.T
    print(f"  cut {b}: covariance rank {int(keep.sum())} of {n}", flush=True)
    sig_out, sig_sc, mean_t, mean_s = [], [], 0.0, 0.0
    S_t = np.zeros(n); S2_t = np.zeros(n); S_s = np.zeros(n); S2_s = np.zeros(n); cnt = 0
    B = 256
    for _ in range(N // (2 * B)):
        w = rng.standard_normal((int(keep.sum()), B))
        for sgn in (1.0, -1.0):                    # antithetic: the odd kappa3 part changes sign, the even kappa4 part not
            ws = sgn * w
            z = mu[:, None] + Lc @ ws
            v = Linv.T @ ws
            U = Pi @ ((k3[:, None] * (v * v - Pd[:, None])) / 6.0 + (k4[:, None] * (v ** 3 - 3.0 * Pd[:, None] * v)) / 24.0)
            h = (np.sum(k3[:, None] * (v ** 3 - 3.0 * Pd[:, None] * v), axis=0) / 6.0
                 + np.sum(k4[:, None] * (v ** 4 - 6.0 * Pd[:, None] * v * v + 3.0 * Pd[:, None] ** 2), axis=0) / 24.0)
            # suffix: y_b = relu(z_b), then layers b+1..L-1; tangent t through the actual gates
            g = (z > 0); y = z * g; t = U * g
            for l in range(b + 1, L):
                zl = Wcol[l] @ y; tl = Wcol[l] @ t
                g = (zl > 0); y = zl * g; t = tl * g
            S_t += t.sum(1); S2_t += (t * t).sum(1)
            hy = y * h[None, :]
            S_s += hy.sum(1); S2_s += (hy * hy).sum(1); cnt += B
    mt_ = S_t / cnt; vt = S2_t / cnt - mt_ ** 2
    ms_ = S_s / cnt; vs = S2_s / cnt - ms_ ** 2
    rms = lambda x: float(np.sqrt(np.mean(x)))
    sig = rms(mt_ ** 2); sdt = rms(vt); sds = rms(vs)
    need = (sdt / (0.1 * sig)) ** 2 if sig > 0 else float("nan")
    fl = need * (L - b) * 2 * 2 * n * n           # value + tangent passes through the suffix
    print(f"net {net} cut {b:2d} (N {cnt}): signal rms {sig:.3e} | per-sample sd: transport {sdt:.3e}, score {sds:.3e} "
          f"(variance ratio {sds ** 2 / max(sdt ** 2, 1e-300):.0f}) | score-transport mean agreement "
          f"{rms((ms_ - mt_) ** 2) / sig:.2f} of the signal (score noise {sds / math.sqrt(cnt) / sig:.2f}) | samples for 10%: "
          f"{need:.0f}, suffix FLOPs {fl:.2e} ({fl / 2 ** 41:.3f} B) | {time.time() - t0:.0f}s", flush=True)
