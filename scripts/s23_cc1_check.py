"""CC1 against the Gaussian closure and Monte Carlo (python scripts/s23_cc1_check.py S21C DATA NET ...).
Per layer: collective moments E c, Var c, Cov(c, r^2), E r^2 (CC1, closure, Monte Carlo) and the mean error against
truth; final layer: RMS error and its coherent part (regression on a_k = w_k . mu_hat)."""
import sys, os, glob, time, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from whest.cc1 import CC1
from scripts.s21_collective import MIX
rms = lambda x: np.sqrt(np.mean(np.square(x)))
SC, D = sys.argv[1], sys.argv[2]
for net in map(int, sys.argv[3:]):
    fs = sorted(glob.glob(f"{SC}/s21c_{net}_*.npz")); NN = 0; mix = 0
    for f in fs: z = np.load(f); NN += int(z["N"]); mix = mix + z["mix"]
    M = {ij: mix[:, kk] / NN for kk, ij in enumerate(MIX)}
    mc = dict(cbar=M[(1, 0)], Vc=M[(2, 0)] - M[(1, 0)] ** 2, cov=M[(1, 1)] - M[(1, 0)] * M[(0, 1)], Er2=M[(0, 1)])
    W = [np.asarray(w, dtype=np.float64) for w in np.load(f"{D}/W_off{net}.npy")]
    T = np.load(f"{D}/truth_off{net}.npz")["m"].astype(np.float64)
    a = W[-1] @ (T[-2] / np.linalg.norm(T[-2])); x = a / a.std()
    res = {}
    for name, model in (("closure", CC1(collective=False)), ("cc1", CC1(Q=7))):
        t0 = time.time(); out, coll = model.run(W); res[name] = (out, coll, time.time() - t0)
    print(f"\n=== net {net}")
    print(" l | E c rel err (closure, cc1) | Var c / MC - 1 (closure, cc1) | Cov(c,r^2) / MC - 1 (closure, cc1) | E r^2 rel err | mean RMS err (closure, cc1)")
    for l in (0, 1, 2, 4, 8, 12, 14, 15):
        cl, cc = res["closure"][1][l], res["cc1"][1][l]
        f = lambda v, k: v / mc[k][l] - 1
        print(f"{l:2d} | {f(cl[0], 'cbar'):+.1e} {f(cc[0], 'cbar'):+.1e} | {f(cl[1], 'Vc'):+.3f} {f(cc[1], 'Vc'):+.3f} | {f(cl[2], 'cov'):+.3f} {f(cc[2], 'cov'):+.3f} |"
              f" {f(cl[3], 'Er2'):+.1e} {f(cc[3], 'Er2'):+.1e} | {rms(res['closure'][0][l] - T[l]):.2e} {rms(res['cc1'][0][l] - T[l]):.2e}")
    for name in ("closure", "cc1"):
        e = res[name][0][-1] - T[-1]; coh = np.polyval(np.polyfit(x, e, 8), x)
        print(f"  {name:8s}: final RMS {rms(e):.3e}, coherent {rms(coh):.2e}, incoherent {rms(e - coh):.2e}  ({res[name][2]:.0f}s)", flush=True)
