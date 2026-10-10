# Note XLIV section 9: the chain's own effective gain charges (per layer, per slice), built exactly as the V62_GPK hook builds them
# (forms from the chain's own mean, variance and off-diagonal covariance), from a baseline chain dump.
#   python -I pin_targets.py LOCDIR CD2DIR NET OUT.npz
import sys, numpy as np
loc, cdd, net, outp = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
W = np.load(f"{loc}/W_off{net}.npy"); c2 = np.load(f"{cdd}/chaindump2_{net}.npz"); out = c2["out"].astype(np.float64)
n = W.shape[1]; off = ~np.eye(n, dtype=bool); f64 = lambda a: np.asarray(a, dtype=np.float64)
T = {k: np.full(16, np.nan) for k in ("e_v", "e_3", "e_21", "e_22", "e_31", "e_4")}
for s in range(1, 15):
    mu = f64(W[s]) @ out[s - 1]; v = f64(c2[f"var_{s}"]); Co = f64(c2[f"C_off_{s}"]); np.fill_diagonal(Co, 0.0)
    T["e_v"][s] = float(mu @ (Co @ mu + v * mu) / (mu @ mu) ** 2)
    f3 = 6 * mu * v; f4 = 12 * v * v
    f21 = np.where(off, 2 * (2 * mu[:, None] * Co + mu[None, :] * v[:, None]), 0.0)
    f22 = np.where(off, 4 * v[:, None] * v[None, :] + 8 * Co * Co, 0.0)
    f31w = np.where(off, 12 * Co * v[None, :], 0.0)
    T["e_3"][s] = float(f3 @ f64(c2[f"D3_{s}"]) / (f3 @ f3)); T["e_4"][s] = float(f4 @ f64(c2[f"g4row_{s}"]) / (f4 @ f4))
    T["e_21"][s] = float(np.sum(f21 * f64(c2[f"D21_{s}"])) / np.sum(f21 * f21))
    T["e_22"][s] = float(np.sum(f22 * f64(c2[f"wk4m_{s}"])) / np.sum(f22 * f22))
    T["e_31"][s] = float(np.sum(f31w * f64(c2[f"wk431_{s}"])) / np.sum(f31w * f31w))
np.savez(outp, **T)
for s in (3, 8, 14):
    print(f"net {net} layer {s}: " + " ".join(f"{k} {1e3 * T[k][s]:6.2f}" for k in T))
