"""Gaussian closure vs the Monte Carlo cache (S = 2^20): conventions of the cache: s1..s4 = raw moments of the
PRE-activation z_{l+1}; Hc = E[z z^T]; pos = P(z > 0); xs1 = sum of post-activations; G = E[x h_{l+1}^T] (n_in x n).
  python scripts/diag_gauss_mc.py DATA_DIR NET
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from whest.relu_gauss import gauss_chain
from scipy.special import ndtr
D, net = sys.argv[1], int(sys.argv[2])
W = np.load(f"{D}/W_off{net}.npy").astype(np.float64); mc = np.load(f"{D}/mccache_{net}.npz"); S = float(mc["S"])
mt = np.load(f"{D}/truth_off{net}.npz")["m"].astype(np.float64)
rec = {}; out = gauss_chain(W, K=10, record=rec); L, n, _ = W.shape
z1 = mc["s1"] / S; z2 = mc["s2"] / S; z3 = mc["s3"] / S; var_mc = z2 - z1 ** 2
k3_mc = z3 - 3 * z2 * z1 + 2 * z1 ** 3
hpost = mc["xs1"] / S
Lprev = np.eye(n)
print("layer | pre-mean rms err (floor) | pre-var relerr | off-diag cov relF | post-mean rms err: chain vs MC, MC vs truth, chain vs truth | 1st-chaos share post (med/mean) | Lmg 1-step relF (noise est) | skew(z) rms | gate err rms")
for l in range(L):
    r = rec[l]; Cz = mc["Hc"][l].astype(np.float64) - np.outer(z1[l], z1[l])
    me = np.sqrt(np.mean((r["mu"] - z1[l]) ** 2)); floor = np.sqrt(np.mean(var_mc[l]) / S)
    dv = np.sqrt(np.mean(((r["sigma"] ** 2 - var_mc[l]) / var_mc[l]) ** 2))
    offm = Cz - np.diag(np.diag(Cz)); offc = r["C"] - np.diag(np.diag(r["C"]))
    offe = np.linalg.norm(offc - offm) / np.linalg.norm(offm)
    pm_cm = np.sqrt(np.mean((r["m"] - hpost[l]) ** 2)); pm_mt = np.sqrt(np.mean((hpost[l] - mt[l]) ** 2)); pm_ct = np.sqrt(np.mean((r["m"] - mt[l]) ** 2))
    G = mc["G"][l].astype(np.float64)
    varpost = np.diag(r["Kh"])     # chain's post variance (close enough for a share)
    share = np.sum(G ** 2, 0) / varpost
    Lh = (Lprev @ W[l].T) * ndtr(r["alpha"])[None, :]
    lerr = np.linalg.norm(Lh - G) / np.linalg.norm(G)
    noise = np.sqrt(n * n * np.mean(r["Kh"].diagonal() + r["m"] ** 2) / S) / np.linalg.norm(G)
    skew = np.sqrt(np.mean((k3_mc[l] / var_mc[l] ** 1.5) ** 2))
    gate = np.sqrt(np.mean((ndtr(r["alpha"]) - mc["pos"][l] / S) ** 2))
    print(f"{l:5d} | {me:8.1e} ({floor:7.1e}) | {dv:8.1e} | {offe:8.3f} | {pm_cm:8.1e} {pm_mt:8.1e} {pm_ct:8.1e} | {np.median(share):6.3f}/{np.mean(share):6.3f} | {lerr:6.4f} ({noise:6.4f}) | {skew:7.4f} | {gate:7.1e}")
    Lprev = G
