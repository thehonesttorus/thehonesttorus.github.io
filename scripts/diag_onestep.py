"""One-step closure tests for the collective (top) mode, with Monte Carlo inputs.
  python scripts/diag_onestep.py DATA NET MCS_ALL.npz [opt=val]
Uses: mccache (mu, E[zz^T] for every layer, S=2^20) for targets; MCS (16M samples) for the slice inputs D21, K22, K31
at its layers; the chain (v3) for comparison of its own D21/K22 against MC. For each MCS layer l (pre-activation z_l):
  chain-vs-MC: D21 corr/slope, K22 corr/slope and its gain projection 4v = s^T K22 s / (s^T s)^2 (s = E z^2);
  one-step: from MC (mu_l, C_l[, D21, K22, K31]) build Cov(h_l) by the chain's closure and compare W C W^T with MC C_{l+1}:
  top-eigenvalue ratio, mu-direction ratio, trace ratio, for the Gaussian closure, +kappa3 (D21), +kappa4 (2,2), +(3,1), +kappa3/4 in E[h^2].
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from scipy.special import ndtr
from whest.k3chain3 import k3_chain3
from whest.relu_gauss import phi, hermite_relu
D, net, mcsf = sys.argv[1], int(sys.argv[2]), sys.argv[3]; opts = {}
for a in sys.argv[4:]:
    k, v = a.split("="); opts[k] = v if k in ("k4", "k22gate") else (float(v) if "." in v else int(v))
W = np.load(f"{D}/W_off{net}.npy").astype(np.float64); mc = np.load(f"{D}/mccache_{net}.npz"); S = float(mc["S"])
ms = np.load(mcsf); N = float(ms["N"]); layers = [int(x) for x in ms["layers"]]
z1 = mc["s1"] / S
rec = {}; out, _ = k3_chain3(W, dict(opts, k22=1), record=rec); L, n, _ = W.shape
def central_slices(l):
    m1 = ms["s1"][l] / N; m2 = ms["s2"][l] / N; M2 = ms[f"M2_{l}"] / N; M21 = ms[f"M21_{l}"] / N; M22 = ms[f"M22_{l}"] / N; M31 = ms[f"M31_{l}"] / N
    C = M2 - np.outer(m1, m1); var = np.diag(C)
    # D21[i,j] = E[(zi-mi)^2 (zj-mj)] = M21 - 2 mi M2 - mj m2_i + 2 mi^2 mj
    D21 = M21 - 2 * m1[:, None] * M2 - m1[None, :] * m2[:, None] + 2 * (m1 ** 2)[:, None] * m1[None, :]
    M12 = M21.T
    cm22 = (M22 - 2 * m1[None, :] * M21 - 2 * m1[:, None] * M12 + (m1 ** 2)[None, :] * m2[:, None] + (m1 ** 2)[:, None] * m2[None, :]
            + 4 * np.outer(m1, m1) * M2 - 2 * np.outer(m1, m1 ** 2) * m1[:, None] - 2 * np.outer(m1 ** 2, m1) * m1[None, :] + np.outer(m1 ** 2, m1 ** 2))
    K22 = cm22 - np.outer(var, var) - 2 * C ** 2
    # E[(zi-mi)^3 (zj-mj)] = M31 - 3 mi M21 + 3 mi^2 M2 - mi^3 m1_j - mj (m3_i - 3 mi m2_i + 2 mi^3) ... build via central moments:
    m3 = ms["s3"][l] / N
    cm31 = (M31 - m1[None, :] * m3[:, None] - 3 * m1[:, None] * M21 + 3 * np.outer(m1, m1) * m2[:, None] + 3 * (m1 ** 2)[:, None] * M2
            - 3 * np.outer(m1 ** 2, m1) * m1[:, None] - (m1 ** 3)[:, None] * m1[None, :] + np.outer(m1 ** 3, m1))
    K31 = cm31 - 3 * var[:, None] * C
    return m1, C, D21, K22, K31
def relu_cov(mu, C, D21=None, K22=None, K31=None, K=8, d3=None, k4=None):
    var = np.clip(np.diag(C), 1e-30, None); sigma = np.sqrt(var); alpha = mu / sigma; Phi = ndtr(alpha); ph = phi(alpha)
    w2 = ph / sigma; w3 = -alpha * ph / var; d = hermite_relu(alpha, K); m = sigma * d[0]
    rho = C / np.outer(sigma, sigma); np.fill_diagonal(rho, 0.0)
    Kh = np.zeros_like(C); rk = np.ones_like(rho); fact = 1.0
    for k in range(1, K + 1):
        rk = rk * rho; fact *= k; Kh += np.outer(d[k], d[k]) * rk / fact
    Kh *= np.outer(sigma, sigma)
    if D21 is not None: Kh += 0.5 * (D21 * np.outer(w2, Phi) + D21.T * np.outer(Phi, w2))
    if K22 is not None: Kh += 0.25 * K22 * np.outer(w2, w2)
    if K31 is not None: Kh += (K31 * np.outer(w3, Phi) + K31.T * np.outer(Phi, w3)) / 6.0
    second = var * ((1 + alpha ** 2) * Phi + alpha * ph)
    if d3 is not None: second = second + d3 * ph / (3 * sigma); m = m - d3 * alpha * ph / (6 * var)
    if k4 is not None: second = second - k4 * alpha * ph / (12 * var); m = m + k4 * (alpha ** 2 - 1) * ph / (24 * sigma ** 3)
    np.fill_diagonal(Kh, second - m * m)
    return Kh, m
print("MCS layers", layers)
for l in layers:
    if l + 1 >= L: continue
    m1, C, D21, K22, K31 = central_slices(l)
    var = np.diag(C); s = var + m1 ** 2
    r = rec[l]; Cm_next = mc["Hc"][l + 1].astype(np.float64) - np.outer(z1[l + 1], z1[l + 1])
    w, Q = np.linalg.eigh(Cm_next); q = Q[:, np.argmax(w)]; lam = w.max(); u = z1[l + 1] / np.linalg.norm(z1[l + 1])
    off = ~np.eye(n, dtype=bool)
    print(f"--- layer {l} (pre-activation z_{l}) ---")
    print(f"  chain vs MC: D21 corr {np.corrcoef(r['D21'][off], D21[off])[0,1]:+.3f} slope {np.dot(r['D21'][off], D21[off])/np.dot(D21[off], D21[off]):+.3f} | "
          f"K22 off corr {np.corrcoef(r['K22'][off], K22[off])[0,1]:+.3f} slope {np.dot(r['K22'][off], K22[off])/np.dot(K22[off], K22[off]):+.3f} | "
          f"gain 4v: MC {s@K22@s/(s@s)**2:+.3e} chain {s@r['K22']@s/(s@s)**2:+.3e} | K22 diag: MC rms {np.sqrt(np.mean(np.diag(K22)**2)):.2e} chain {np.sqrt(np.mean(np.diag(r['K22'])**2)):.2e}")
    k3d = np.diag(D21); k4d = np.diag(K22)
    for name, kw in [("Gaussian", {}), ("+D21", dict(D21=D21)), ("+D21+K22", dict(D21=D21, K22=K22)), ("+D21+K22+K31", dict(D21=D21, K22=K22, K31=K31)),
                     ("+D21+K22+K31+d3,k4 in E[h^2]", dict(D21=D21, K22=K22, K31=K31, d3=k3d, k4=k4d)), ("+D21 + chain K22", dict(D21=D21, K22=r["K22"]))]:
        Kh, m = relu_cov(m1, C, **kw); Cp = W[l + 1] @ Kh @ W[l + 1].T
        print(f"  one-step {name:32s}: top-mode {q@Cp@q/lam-1:+.4f}  mu-dir {u@Cp@u/(u@Cm_next@u)-1:+.4f}  trace {np.trace(Cp)/np.trace(Cm_next)-1:+.4f}  "
              f"mean rms err vs truth-next {np.sqrt(np.mean((W[l+1]@m - (W[l+1]@(mc['xs1'][l]/S)))**2)):.1e}")
