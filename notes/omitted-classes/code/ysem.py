# Semantics of the chain's post-activation arrays (one-step dump, true z inputs) against Monte Carlo y statistics.
#   python ysem.py NET DUMPFILE MCPREFIX
import sys, numpy as np
net, dmp, pre = int(sys.argv[1]), sys.argv[2], sys.argv[3]
ch = np.load(dmp); F = np.load(f"{pre}_full.npz"); n = F["mu"].shape[1]; off = ~np.eye(n, dtype=bool)
def rc(p, t):
    return f"{np.linalg.norm(p - t) / np.linalg.norm(t):.3f}/{np.corrcoef(p.ravel(), t.ravel())[0, 1]:+.3f}"
for l in range(1, 15):
    if f"K21_{l}" not in ch.files:
        continue
    K21 = ch[f"K21_{l}"].astype(np.float64); K3v = ch[f"K3v_{l}"].astype(np.float64); K11 = ch[f"K11_{l}"].astype(np.float64)
    pk1 = ch[f"pk1v_{l}"].astype(np.float64); K2v = ch[f"K2v_{l}"].astype(np.float64)
    Dy = F["D21_y"][l].astype(np.float64); Cy = F["cov_y"][l].astype(np.float64)
    print(f"layer {l:2d}: mean {rc(pk1, F['mu_y'][l])} var {rc(K2v, F['var_y'][l])} k3 {rc(K3v, F['k3_y'][l])} | "
          f"K21 vs D21_y {rc(K21[off], Dy[off])}, vs D21_y^T {rc(K21[off], Dy.T[off])} | K11+K11^T vs cov_y {rc((K11 + K11.T)[off], Cy[off])}, "
          f"K11 vs cov_y {rc(K11[off], Cy[off])}", flush=True)
