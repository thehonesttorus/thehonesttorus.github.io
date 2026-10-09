"""Per-layer injected mean error of a chain run: e_l = out_l - truth_l; injected_l = e_l - Phi_l (W_l e_{l-1}).
  python scripts/injected.py DATA NET OUTFILE.npz   (OUTFILE from run_k3v2.py with OUT set) or pass opts to run inline:
  python scripts/injected.py DATA NET -- k22born=1 ..."""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from scipy.special import ndtr
from whest.k3chain3 import k3_chain3 as k3_chain2
D, net = sys.argv[1], int(sys.argv[2])
W = np.load(f"{D}/W_off{net}.npy").astype(np.float64); mt = np.load(f"{D}/truth_off{net}.npz")["m"].astype(np.float64)
if sys.argv[3] == "--":
    opts = {}
    for a in sys.argv[4:]:
        k, v = a.split("="); opts[k] = v if k in ("k4", "k22gate") else (float(v) if "." in v else int(v))
    rec = {}; out, fl = k3_chain2(W, opts, record=rec); alphas = [rec[l]["alpha"] for l in range(len(rec))]
else:
    z = np.load(sys.argv[3]); out = z["out"]; alphas = None; opts = sys.argv[3]
e = out - mt; L = len(out)
print(f"{opts}: final MSE {np.mean(e[-1]**2):.3e}")
print("layer | rms err | rms injected | injected/total | coherent share of injected (proj on truth-mean direction)")
for l in range(L):
    if l == 0:
        inj = e[0]
    else:
        g = ndtr(alphas[l]) if alphas is not None else 0.5
        inj = e[l] - g * (W[l] @ e[l - 1])
    coh = (np.dot(inj, mt[l]) / np.linalg.norm(mt[l])) ** 2 / np.sum(inj ** 2)
    print(f"{l:5d} | {np.sqrt(np.mean(e[l]**2)):7.1e} | {np.sqrt(np.mean(inj**2)):7.1e} | {np.sum(inj**2)/np.sum(e[l]**2):6.3f} | {coh:6.3f}")
