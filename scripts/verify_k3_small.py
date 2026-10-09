"""Verify the CP third-cumulant sources against brute-force Monte Carlo on a tiny He-init ReLU network.
  python scripts/verify_k3_small.py [n=48] [L=3] [N=4000000] [seed=0]
Reports, per layer: relative rms error of the chain's D3 (diag kappa3 of z) and D21 (kappa3(z_i,z_i,z_j)) against
MC, the MC noise level, and the mean readout error with/without the Edgeworth term.
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from whest.k3chain import k3_chain
from whest.relu_gauss import gauss_chain
kw = dict(n=48, L=3, N=4_000_000, seed=0)
for a in sys.argv[1:]:
    k, v = a.split("="); kw[k] = int(v)
n, L, N, seed = kw["n"], kw["L"], kw["N"], kw["seed"]
rng = np.random.default_rng(seed)
W = (rng.standard_normal((L, n, n)) * np.sqrt(2.0 / n)).astype(np.float64)
# Monte Carlo: accumulate E[z], E[z z^T], E[z_i z_i z_j] (full (2,1) slices), E[h] per layer
B = 100_000; acc = [dict(s1=np.zeros(n), s2=np.zeros((n, n)), s21=np.zeros((n, n)), h=np.zeros(n)) for _ in range(L)]
for _ in range(N // B):
    h = rng.standard_normal((B, n))
    for l in range(L):
        z = h @ W[l].T; h = np.maximum(z, 0)
        a = acc[l]; a["s1"] += z.sum(0); a["s2"] += z.T @ z; a["s21"] += (z * z).T @ z; a["h"] += h.sum(0)
rec = {}; out, _ = k3_chain(W, dict(hub=1, K=8), record=rec); outg = gauss_chain(W, K=8)
for l in range(L):
    a = acc[l]; m1 = a["s1"] / N; m2 = a["s2"] / N; m21 = a["s21"] / N; hm = a["h"] / N
    C = m2 - np.outer(m1, m1)
    # kappa3(z_i,z_i,z_j) = E[z_i^2 z_j] - 2 m_i C_ij - m_j C_ii - m_i^2 m_j ... compute via central moments:
    # E[(zi-mi)^2 (zj-mj)] = m21_ij - 2 m_i m2_ij - m_j m2_ii + 2 m_i^2 m_j
    k21 = m21 - 2 * m1[:, None] * m2 - m1[None, :] * np.diag(m2)[:, None] + 2 * (m1 ** 2)[:, None] * m1[None, :]
    k3 = np.diag(k21).copy()
    r = rec[l]
    D3 = r["D3"]; var = r["var"]
    noise3 = np.sqrt(15 * np.mean(var ** 3) / N)      # std of a sample third central moment ~ sqrt(15 sigma^6 / N) at kurtosis 3
    e3 = np.sqrt(np.mean((D3 - k3) ** 2)); s3 = np.sqrt(np.mean(k3 ** 2))
    print(f"layer {l}: |k3_mc| rms {s3:.3e}  chain D3 err rms {e3:.3e}  (MC noise ~{noise3:.1e})  corr {np.corrcoef(D3, k3)[0,1]:.4f}  "
          f"slope {np.dot(D3, k3)/np.dot(k3, k3):.3f} | mean err: gauss {np.sqrt(np.mean((outg[l]-hm)**2)):.2e}  k3chain {np.sqrt(np.mean((out[l]-hm)**2)):.2e}  (MC floor {np.sqrt(np.mean(np.maximum(hm,0)*0+0.3)/N):.1e})")
    if l < L - 1 and "D21" in r:
        D21 = r["D21"]; off = ~np.eye(n, dtype=bool)
        e21 = np.sqrt(np.mean((D21 - k21)[off] ** 2)); s21 = np.sqrt(np.mean(k21[off] ** 2))
        print(f"         D21 off-diag: |k21_mc| rms {s21:.3e}  err rms {e21:.3e}  corr {np.corrcoef(D21[off], k21[off])[0,1]:.4f}  slope {np.dot(D21[off], k21[off])/np.dot(k21[off], k21[off]):.3f}")
