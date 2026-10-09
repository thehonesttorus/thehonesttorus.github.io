"""Verify chain v2 on a tiny network against Monte Carlo (same protocol as verify_k3_small.py)."""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from whest.k3chain2 import k3_chain2
kw = dict(n=256, L=3, N=3_000_000, seed=0); opts = {}
for a in sys.argv[1:]:
    k, v = a.split("=")
    if k in kw: kw[k] = int(v)
    else: opts[k] = v if k == "k4" else (float(v) if "." in v else int(v))
n, L, N, seed = kw["n"], kw["L"], kw["N"], kw["seed"]
rng = np.random.default_rng(seed); W = (rng.standard_normal((L, n, n)) * np.sqrt(2.0 / n))
B = 100_000; acc = [dict(s1=np.zeros(n), s2=np.zeros((n, n)), s21=np.zeros((n, n)), h=np.zeros(n)) for _ in range(L)]
for _ in range(N // B):
    h = rng.standard_normal((B, n))
    for l in range(L):
        z = h @ W[l].T; h = np.maximum(z, 0)
        a = acc[l]; a["s1"] += z.sum(0); a["s2"] += z.T @ z; a["s21"] += (z * z).T @ z; a["h"] += h.sum(0)
rec = {}; out, _ = k3_chain2(W, opts, record=rec)
print("opts", opts)
for l in range(L):
    a = acc[l]; m1 = a["s1"] / N; m2 = a["s2"] / N; m21 = a["s21"] / N; hm = a["h"] / N
    k21 = m21 - 2 * m1[:, None] * m2 - m1[None, :] * np.diag(m2)[:, None] + 2 * (m1 ** 2)[:, None] * m1[None, :]
    k3 = np.diag(k21).copy(); r = rec[l]; D3 = r["D3"]; var = r["var"]
    noise3 = np.sqrt(6 * np.mean(var ** 3) / N)
    if l == 0:
        print(f"layer 0: Gaussian; mean err {np.sqrt(np.mean((out[0]-hm)**2)):.2e} (MC floor {np.sqrt(0.34*np.mean(var)/N):.1e})"); continue
    print(f"layer {l}: |k3_mc| rms {np.sqrt(np.mean(k3**2)):.3e}  D3 err rms {np.sqrt(np.mean((D3-k3)**2)):.3e} (noise {noise3:.1e})  corr {np.corrcoef(D3,k3)[0,1]:.4f}  slope {np.dot(D3,k3)/np.dot(k3,k3):.3f} | mean err {np.sqrt(np.mean((out[l]-hm)**2)):.2e}")
    if l < L - 1:
        D21 = r["D21"]; off = ~np.eye(n, dtype=bool); noise21 = np.sqrt(2 * np.mean(var) ** 1.5 / N)
        print(f"         D21 off: |k21| rms {np.sqrt(np.mean(k21[off]**2)):.3e}  err rms {np.sqrt(np.mean((D21-k21)[off]**2)):.3e} (noise {noise21:.1e})  corr {np.corrcoef(D21[off],k21[off])[0,1]:.4f}  slope {np.dot(D21[off],k21[off])/np.dot(k21[off],k21[off]):.3f}")
