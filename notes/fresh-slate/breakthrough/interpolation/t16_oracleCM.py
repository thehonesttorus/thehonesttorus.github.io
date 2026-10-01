"""T16: value of (i) the exact post-activation covariance and (ii) a kappa4 matrix-trace channel M inside TC (oracle,
MC from t15), MLP 0 at n=1024.  Chaos parts of kappa4(z) = kappa4(a)[w...]:
  kappa4(z_p)        = 3 s4 X + 6 s2 mh_p,         mh_p = (W^T M W)_pp - s2 tr M
  kappa4(z_p,z_p,z_q,z_q) = s4 X + s2 (mh_p + mh_q)  (p != q)
  kappa4(z_p,z_p,z_p,z_q) = 3 s2 (W^T M W)_pq
Second-order Edgeworth: dE a = k4 phi (a^2-1)/(24 v s), dE a^2 = -k4 a phi/(12 v),
dCov_pq = 1/4 k4(ppqq) (phi/s)_p (phi/s)_q + 1/6 k4(pppq) (-a phi/v)_p Phi_q + (p<->q)."""
import sys, numpy as np
from common import bench, phi, Phi, chi_mean_ratio
from t2_closures import hk
import tc
f = sys.argv[1]
S_ = bench.load_set("w1024_d16"); W = bench.weights(S_, 0).astype(np.float64); T = S_["means"][0]; nz = S_["noise"][0]
D = np.load(f); Cmc = D["C"].astype(np.float64); Mmc = D["M"].astype(np.float64); Xmc = D["X"]
n = 1024; L = 16; s2 = 2.0 / n; s4 = s2 * s2
def run(trueC=False, Mch=False, Xonly=False, k4scale=1.0):
    _, ts = tc.predict(W, ret_state=True)
    out = []; mu = None; C = None
    for l in range(L):
        Wl = W[l]
        if l == 0: m = np.zeros(n); S = Wl.T @ Wl; y = np.zeros(n)
        else:
            m = mu @ Wl; S = Wl.T @ (Cmc[l - 1] if trueC else C) @ Wl; y = Wl.T @ ts[l - 1]
        v = np.diag(S); s = np.sqrt(v); a = m / s; P = Phi(a); p = phi(a)
        k3 = 3 * s2 * y
        mu_n = m * P + s * p - k3 * a * p / (6 * v); sec = (m*m+v)*P + m*s*p + k3 * p / (3 * s)
        R = S / np.outer(s, s); H = hk(a, 2)
        Cn = np.outer(s*H[0], s*H[0]) * R + 0.5*np.outer(s*H[1], s*H[1]) * R * R
        u1 = p / s; Cn += 0.5 * s2 * (np.outer(u1, P * y) + np.outer(P * y, u1))
        if (Mch or Xonly) and l > 0:
            Ml = Mmc[l - 1]; X = Xmc[l - 1]
            if Xonly: WMW = np.zeros((n, n)); mh = np.zeros(n)
            else: WMW = Wl.T @ Ml @ Wl; mh = np.diag(WMW) - s2 * np.trace(Ml)
            k4d = k4scale * (3 * s4 * X + 6 * s2 * mh)
            mu_n = mu_n + k4d * p * (a * a - 1) / (24 * v * s); sec = sec - k4d * a * p / (12 * v)
            k22 = k4scale * (s4 * X + s2 * (mh[:, None] + mh[None, :]))
            k31 = k4scale * 3 * s2 * WMW
            e3 = -a * p / v
            Cn += 0.25 * k22 * np.outer(u1, u1) + (k31 * e3[:, None] * P[None, :] + (k31 * e3[:, None] * P[None, :]).T) / 6
        np.fill_diagonal(Cn, np.maximum(sec - mu_n ** 2, 1e-12)); mu, C = mu_n, Cn; out.append(mu)
    return np.stack(out) * chi_mean_ratio(n)
for kw in [dict(), dict(trueC=True), dict(Xonly=True), dict(Mch=True), dict(Mch=True, k4scale=0.5)]:
    p = run(**kw); print(kw, f"final raw {((p[-1]-T[-1])**2).mean()-nz:.3e}  bias {np.mean(p[-1]-T[-1]):+.1e}  per-layer:",
                         " ".join(f"{x:.1e}" for x in ((p-T)**2).mean(1)[[3,7,11,15]]), flush=True)
