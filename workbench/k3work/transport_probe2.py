# The transport estimator for a physical old kappa3 source, against a benign control (note ray-compiler, section 5).
#   python transport_probe2.py NET MCPREFIX CUT AGE [--n N]
# Source "legs": the diagonal class of the source born at the post-activation of layer b0 = CUT - AGE, carried to the
# pre-activation of CUT by the first-order gated transport, T = sum_i s_i l_i x l_i x l_i, s = kappa3(y_b0) diagonal
# (Monte Carlo), l_i = columns of L = W_CUT diag(Phi_(CUT-1)) W_(CUT-1) ... diag(Phi_b0) (true gates Phi = Phi(mu/sigma)).
# Its Stein field at the reference N(mu, Sigma) of z_CUT (support-restricted) is U = (1/6) L (s o ((L^T v)^2 - diag(L^T P L))),
# v = P (z - mu), and its first-order effect on the output means is E[U . grad F] = E[h F], h = (1/6) sum s_i He3(l_i).
# Control "white": a source diagonal in the reference's whitened coordinates, T = sum_k t (Lc e_k)^(x3), the benign case.
# Reports signal rms (transport mean), per-sample sd of both estimators, their mean agreement, the samples for 10% at the
# rms level, and also |contribution of the source to D3 at CUT| against the true D3 there (a size check).
import sys, math, time, numpy as np
from math import erf
net, pre, b, age = int(sys.argv[1]), sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
N = int(sys.argv[sys.argv.index("--n") + 1]) if "--n" in sys.argv else 4096
Wcol = np.load(f"../official/W_off{net}.npy").astype(np.float64)
Ln, n, _ = Wcol.shape
F = np.load(f"{pre}_full.npz")
mu_all, var_all = F["mu"].astype(np.float64), F["var"].astype(np.float64)
Phi = 0.5 * (1 + np.vectorize(erf)(mu_all / np.sqrt(var_all) / math.sqrt(2)))
b0 = b - age
Lg = np.eye(n)
for l in range(b0, b):                              # L = W_b diag(Phi_(b-1)) ... W_(b0+1) diag(Phi_b0)
    Lg = Wcol[l + 1] @ (Phi[l][:, None] * Lg)
s = F["k3_y"][b0].astype(np.float64)
mu = mu_all[b]; S = F["cov"][b].astype(np.float64); S = 0.5 * (S + S.T)
lam, Q = np.linalg.eigh(S); keep = lam > 1e-6 * lam.max(); Qr, lr = Q[:, keep], lam[keep]
Lc = Qr * np.sqrt(lr)[None, :]; Linv = (Qr / np.sqrt(lr)[None, :]).T; Pi = Qr @ Qr.T
P = (Qr / lr[None, :]) @ Qr.T
LPL = np.einsum("ai,ab,bi->i", Lg, P, Lg)          # diag(L^T P L)
d3src = (Lg ** 3) @ s
print(f"net {net} cut {b} age {age}: covariance rank {int(keep.sum())}; |source D3 at cut| {np.linalg.norm(d3src):.3e} "
      f"against |true D3| {np.linalg.norm(F['k3'][b]):.3e}; leg conditioning: median l^T P l / |l|^2 "
      f"{np.median(LPL / np.sum(Lg * Lg, axis=0)):.3e} (1/median var {1 / np.median(np.diag(S)):.3e})", flush=True)
# benign control: whitened diagonal with the same typical diagonal size at z
t_white = float(np.median(np.abs(d3src)) / np.median(np.sum(Lc ** 3, axis=1) ** 2) ** 0.5) if False else 0.0
rng = np.random.default_rng(5)
res = {}
for name in ("legs", "white"):
    St = np.zeros(n); S2t = np.zeros(n); Ss = np.zeros(n); S2s = np.zeros(n); cnt = 0
    Sc = np.zeros(n); S2c = np.zeros(n); Stc = np.zeros(n); Sd = np.zeros(n); S2d = np.zeros(n)
    for _ in range(N // 512):
        w = rng.standard_normal((Lc.shape[1], 256))
        for sg in (1.0, -1.0):
            ws = sg * w; z = mu[:, None] + Lc @ ws; v = Linv.T @ ws
            if name == "legs":
                u = Lg.T @ v                                           # (n, B)
                U = Pi @ (Lg @ (s[:, None] * (u * u - LPL[:, None]))) / 6.0
                h = np.sum(s[:, None] * (u ** 3 - 3.0 * LPL[:, None] * u), axis=0) / 6.0
            else:                                                      # T = sum_k (Lc e_k)^3 / 6: U = Lc (w^2 - 1)/6
                U = Lc @ (ws * ws - 1.0) / 6.0
                h = np.sum(ws ** 3 - 3.0 * ws, axis=0) / 6.0
            g = z > 0; y = z * g; t = U * g; tc = U * Phi[b][:, None]     # tc: the same tangent through the expected gates
            for l in range(b + 1, Ln):
                zl = Wcol[l] @ y; tl = Wcol[l] @ t; g = zl > 0; y = zl * g; t = tl * g
                tc = (Wcol[l] @ tc) * Phi[l][:, None]
            St += t.sum(1); S2t += (t * t).sum(1); hy = y * h[None, :]; Ss += hy.sum(1); S2s += (hy * hy).sum(1); cnt += 256
            Sc += tc.sum(1); S2c += (tc * tc).sum(1); Stc += (t * tc).sum(1)
            d1 = t - tc; Sd += d1.sum(1); S2d += (d1 * d1).sum(1)
    mt_, ms_ = St / cnt, Ss / cnt
    vt, vs = S2t / cnt - mt_ ** 2, S2s / cnt - ms_ ** 2
    rms = lambda x: float(np.sqrt(np.mean(x)))
    sig, sdt, sds = rms(mt_ ** 2), rms(vt), rms(vs)
    noise_t = sdt / math.sqrt(cnt)
    mc_ = Sc / cnt; vc = S2c / cnt - mc_ ** 2; cov = Stc / cnt - mt_ * mc_
    beta = cov / np.maximum(vc, 1e-300); vopt = vt - beta * cov        # optimal per-neuron control coefficient
    vd = S2d / cnt - (Sd / cnt) ** 2                                     # coefficient 1 (no fit)
    print(f"  {name:5s} expected-gate control: E[control] rms {rms(mc_ ** 2):.3e} (exact 0; noise {rms(vc) / math.sqrt(cnt):.3e}); "
          f"per-sample sd with coefficient 1 {rms(vd):.3e}, with optimal coefficients {rms(np.maximum(vopt, 0)):.3e} "
          f"(median beta {np.median(beta):.2f}); variance reduction x{sdt ** 2 / max(rms(vd) ** 2, 1e-300):.1f} / "
          f"x{sdt ** 2 / max(rms(np.maximum(vopt, 0)) ** 2, 1e-300):.1f}", flush=True)
    print(f"  {name:5s} (N {cnt}): transport mean rms {sig:.3e} (its noise {noise_t:.3e}) | per-sample sd: transport {sdt:.3e}, "
          f"score {sds:.3e} (variance ratio {sds ** 2 / max(sdt ** 2, 1e-300):.0f}) | score - transport mean rms {rms((ms_ - mt_) ** 2):.3e} "
          f"(score noise {sds / math.sqrt(cnt):.3e}) | samples for 10% of the mean: {(sdt / (0.1 * max(sig - noise_t, 1e-30))) ** 2:.0f}", flush=True)
