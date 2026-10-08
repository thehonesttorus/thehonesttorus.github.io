# Size of the kappa3-driven terms vs the one-step residual: readout skew term (F3/6) D3 at layer l, and the
# next-layer variance moved by the D21 Wick term of the post covariance at l-1 (projected with the fresh rows).
import pickle, numpy as np
from scipy.stats import norm
D = {d["layer"]: d for d in pickle.load(open("dump_oracleall_off0.pkl", "rb"))}
pred = np.load("pred_oracleall_off0.npy"); mt = np.load("../official/truth_off0.npz")["m"].astype(np.float64)
Wc = np.load("../official/W_off0.npy").astype(np.float64)
for l in sorted(D):
    c = D[l]; s = np.sqrt(c["var"]); a = c["mu"] / s; ph = norm.pdf(a)
    skew = (-a * ph / s**2) / 6 * c["D3"]
    line = f"layer {l:2d}: rms r {np.sqrt(np.mean((mt[l]-pred[l])**2)):.1e} | readout skew term {np.sqrt(np.mean(skew**2)):.1e}"
    if l - 1 in D and D[l - 1]["D21"] is not None:
        p = D[l - 1]; sp = np.sqrt(p["var"]); ap = p["mu"] / sp; F1, F2 = norm.cdf(ap), norm.pdf(ap) / sp
        M = np.outer(F1, F2) * p["D21"].T; M = 0.5 * (M + M.T); np.fill_diagonal(M, 0)
        dv = np.einsum("ij,ij->i", Wc[l] @ M, Wc[l]) * ph / (2 * s)
        line += f" | d21-moved mean {np.sqrt(np.mean(dv**2)):.1e}"
    print(line)
