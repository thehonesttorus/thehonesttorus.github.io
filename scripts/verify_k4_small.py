"""Verify chain v3's fourth-cumulant pieces on a tiny network against Monte Carlo (diagonal kappa4 of z and the (2,2)
slice).   python scripts/verify_k4_small.py n=256 L=3 N=6000000 [opts]"""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from whest.k3chain3 import k3_chain3
kw = dict(n=256, L=3, N=6_000_000, seed=0); opts = {}
for a in sys.argv[1:]:
    k, v = a.split("=")
    if k in kw: kw[k] = int(v)
    else: opts[k] = v if k in ("k4", "k22gate") else (float(v) if "." in v else int(v))
n, L, N, seed = kw["n"], kw["L"], kw["N"], kw["seed"]
rng = np.random.default_rng(seed); W = (rng.standard_normal((L, n, n)) * np.sqrt(2.0 / n))
B = 50_000; acc = [dict(s1=np.zeros(n), s2=np.zeros(n), s3=np.zeros(n), s4=np.zeros(n), M2=np.zeros((n, n)), M22=np.zeros((n, n)), M12=np.zeros((n, n)), h=np.zeros(n)) for _ in range(L)]
for _ in range(N // B):
    h = rng.standard_normal((B, n))
    for l in range(L):
        z = h @ W[l].T; h = np.maximum(z, 0); a = acc[l]; z2 = z * z
        a["s1"] += z.sum(0); a["s2"] += z2.sum(0); a["s3"] += (z2 * z).sum(0); a["s4"] += (z2 * z2).sum(0)
        a["M2"] += z.T @ z; a["M22"] += z2.T @ z2; a["M12"] += z.T @ z2; a["h"] += h.sum(0)
rec = {}; out, _ = k3_chain3(W, dict(opts, k22=1), record=rec)
print("opts", opts)
for l in range(1, L):
    a = acc[l]; m1 = a["s1"] / N; m2 = a["s2"] / N; m3 = a["s3"] / N; m4 = a["s4"] / N; var = m2 - m1 ** 2
    k4 = m4 - 4 * m3 * m1 + 6 * m2 * m1 ** 2 - 3 * m1 ** 4 - 3 * var ** 2
    M2 = a["M2"] / N; M22 = a["M22"] / N; M12 = a["M12"] / N      # E[zi zj], E[zi^2 zj^2], E[zi zj^2]
    Cov = M2 - np.outer(m1, m1)
    # central: E[(zi-ai)^2 (zj-aj)^2] = M22 - 2 aj M21_ij - 2 ai M12_ij + aj^2 m2_i + ai^2 m2_j + 4 ai aj M2_ij - 3 ai^2 aj^2 ... derive:
    # (x-a)^2 (y-b)^2 = x^2y^2 - 2b x^2 y - 2a x y^2 + b^2 x^2 + a^2 y^2 + 4ab xy - 2a b^2 x - 2 a^2 b y + a^2 b^2
    M21 = M12.T   # E[zi^2 zj]
    cm22 = (M22 - 2 * m1[None, :] * M21 - 2 * m1[:, None] * M12 + (m1 ** 2)[None, :] * m2[:, None] + (m1 ** 2)[:, None] * m2[None, :]
            + 4 * np.outer(m1, m1) * M2 - 2 * (m1[:, None] * (m1 ** 2)[None, :]) * m1[:, None] - 2 * ((m1 ** 2)[:, None] * m1[None, :]) * m1[None, :] + np.outer(m1 ** 2, m1 ** 2))
    k22 = cm22 - np.outer(var, var) - 2 * Cov ** 2
    r = rec[l]; off = ~np.eye(n, dtype=bool)
    noise4 = np.sqrt(96 * np.mean(var ** 4) / N)
    print(f"layer {l}: kappa4 diag: MC rms {np.sqrt(np.mean(k4**2)):.3e} (noise {noise4:.1e}) | K4_22 corr {np.corrcoef(r['K4_22'],k4)[0,1]:+.3f} slope {np.dot(r['K4_22'],k4)/np.dot(k4,k4):+.3f} | s211 corr {np.corrcoef(r['K4_s211'],k4)[0,1]:+.3f} slope {np.dot(r['K4_s211'],k4)/np.dot(k4,k4):+.3f} | s1111 corr {np.corrcoef(r['K4_s1111'],k4)[0,1]:+.3f} slope {np.dot(r['K4_s1111'],k4)/np.dot(k4,k4):+.3f} | all three: slope {np.dot(r['K4_22']+r['K4_s211']+r['K4_s1111'],k4)/np.dot(k4,k4):+.3f}")
    K22c = r["K22"]
    print(f"         (2,2) off-diag: MC rms {np.sqrt(np.mean(k22[off]**2)):.3e} (noise ~{np.sqrt(np.mean(var**4)/N):.1e}) | chain corr {np.corrcoef(K22c[off],k22[off])[0,1]:+.3f} slope {np.dot(K22c[off],k22[off])/np.dot(k22[off],k22[off]):+.3f} | top-mode: u=mean dir: MC {m1@k22@m1/np.dot(m1,m1)**2:+.3e} chain {m1@K22c@m1/np.dot(m1,m1)**2:+.3e}")
