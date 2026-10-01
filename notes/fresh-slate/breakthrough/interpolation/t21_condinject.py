"""T21: TC chain (MLP 0, n=1024) with MC t (and t^mu) injected: plain chaos-1 vs Bolthausen-conditioned injection
(diag kappa3 and row-scaled (2,1) slice kappa3(z_a,z_a,z_b) ~ r_a^2 ym_perp_b + s2 (y_perp - ym_perp)_b)."""
import numpy as np
from common import bench, phi, Phi, chi_mean_ratio
from t2_closures import hk
S_ = bench.load_set("w1024_d16"); W = bench.weights(S_, 0).astype(np.float64); T = S_["means"][0]; nz = S_["noise"][0]
d = np.load("results/t17_w1024_0.npz"); n = 1024; s2 = 2.0 / n
def run(cond):
    out = []; mu = None; C = None
    for l in range(16):
        Wl = W[l]
        if l == 0: m = np.zeros(n); S = Wl.T @ Wl
        else: m = mu @ Wl; S = Wl.T @ C @ Wl
        v = np.diag(S); s = np.sqrt(v); a = m / s; P = Phi(a); p = phi(a); g = p / s
        if l == 0: k3 = np.zeros(n); rowc = np.zeros(n); colA = np.zeros(n); colB = np.zeros(n)
        else:
            t, tm, muP = d["t"][l - 1], d["tm"][l - 1], d["mu"][l - 1]; nm = np.linalg.norm(muP); muh = muP / nm
            y = Wl.T @ t
            if cond:
                r = (muP @ Wl) / nm; ym = Wl.T @ tm; K3 = tm @ muh; tmu_ = t @ muh
                yp = y - r * tmu_; ymp = ym - r * K3
                k3 = r ** 3 * K3 + 3 * r ** 2 * ymp + 3 * r * s2 * (tmu_ - K3) + 3 * s2 * (yp - ymp)
                rowc = r ** 2 * CONDCOV; colA = ymp; colB = s2 * (yp - ymp) if CONDCOV else s2 * y
            else:
                k3 = 3 * s2 * y; rowc = np.zeros(n); colA = np.zeros(n); colB = s2 * y
        mu_n = m * P + s * p - k3 * a * p / (6 * v); sec = (m*m+v)*P + m*s*p + k3 * p / (3 * s)
        R = S / np.outer(s, s); H = hk(a, 2)
        Cn = np.outer(s*H[0], s*H[0]) * R + 0.5*np.outer(s*H[1], s*H[1]) * R * R
        # kappa_aab = rowc_a colA_b + colB_b ; dC_ab = 1/2 [kappa_aab g_a Phi_b + kappa_bba g_b Phi_a]
        K = np.outer(g * rowc, P * colA) + np.outer(g, P * colB)
        Cn += 0.5 * (K + K.T)
        np.fill_diagonal(Cn, sec - mu_n ** 2); mu, C = mu_n, Cn; out.append(mu)
    return np.stack(out) * chi_mean_ratio(n)
for c, CONDCOV in ((False, 0), (True, 0), (True, 1)):
    p = run(c); print("conditioned" if c else "plain", "cov-cond" if CONDCOV else "", f"final raw {((p[-1]-T[-1])**2).mean()-nz:.3e}")
