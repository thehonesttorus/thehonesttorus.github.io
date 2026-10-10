import sys, glob, numpy as np
import os; sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from scipy.special import ndtr
R, D = sys.argv[1], sys.argv[2]
rms = lambda x: np.sqrt(np.mean(np.square(x)))
for net in map(int, sys.argv[3:]):
    k = np.load(f"{R}/kikm_{net}.npz"); g = np.load(f"{R}/kcg_{net}.npz")
    t, q, tr = k["t"], k["q"], k["truth"]; y_cl = g["m"][-1]; e = y_cl - tr
    W = np.load(f"{D}/W_off{net}.npy")[-1].astype(np.float64); T = np.load(f"{D}/truth_off{net}.npz")["m"]; mh = T[-2] / np.linalg.norm(T[-2]); a = W @ mh; x = a / a.std()
    sd = np.sqrt(q - t * t); al = t / sd; ph = np.exp(-al * al / 2) / np.sqrt(2 * np.pi)
    basis = {"dt (Phi)": ndtr(al), "dvar (phi/2sd)": ph / (2 * sd), "k3 (sd d3)": sd * (-al * ph), "k4 (sd d4)": sd * (al * al - 1) * ph}
    coh = np.polyval(np.polyfit(x, e, 8), x)
    B = np.stack(list(basis.values()), 1); c, *_ = np.linalg.lstsq(B, coh, rcond=None)
    # per-neuron decomposition with the exact contractions: the closure's t, var errors vs MC, and the Gaussian-shape defect
    tg, qg = g["t"][-1], g["q"][-1]; vg = qg - tg * tg; v = q - t * t
    lin = ndtr(al) * (tg - t) + ph / (2 * sd) * (vg - v)
    print(f"net {net}: coherent err {rms(coh):.2e}; fit on [Phi, phi/2sd, sd d3, sd d4] coeffs " + " ".join(f"{cc:+.2e}" for cc in c) +
          f" resid {rms(coh - B @ c):.2e} | closure dt rms {rms(tg - t):.2e} (coherent {rms(np.polyval(np.polyfit(x, tg - t, 8), x)):.2e}),"
          f" dvar rms {rms(vg - v):.2e} (coherent {rms(np.polyval(np.polyfit(x, vg - v, 8), x)):.2e}); linearised err {rms(lin):.2e}, coherent {rms(np.polyval(np.polyfit(x, lin, 8), x)):.2e}")
