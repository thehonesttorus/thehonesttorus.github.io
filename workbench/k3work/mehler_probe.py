# Size of the Gaussian Mehler orders k>=3 missing from the chain's (1,1) slice, from the chain's own per-layer state.
import pickle, numpy as np, math
from scipy.stats import norm
D = pickle.load(open("dump_off0_win0.pkl", "rb"))
Wcol = np.load("../official/W_off0.npy").astype(np.float64)
print(len(D), list(D[0].keys()))
def He(k, a):
    h0, h1 = np.ones_like(a), a
    if k == 0: return h0
    for j in range(1, k):
        h0, h1 = h1, a*h1 - j*h0
    return h1
for d in D:
    li = d.get("layer", None)
    C = np.asarray(d["C_pre"], dtype=np.float64); mu = np.asarray(d["mu"], dtype=np.float64)
    var = np.asarray(d["var"]); s = np.sqrt(var); a = mu/s
    np.fill_diagonal(C, 0.0)
    Co = C; rho = Co/np.outer(s, s)
    ph = norm.pdf(a)
    out = []
    for k in range(2, 8):
        Fk = He(k-2, a)*ph/s**(k-1)
        T = np.outer(Fk, Fk)*Co**k/math.factorial(k)
        out.append(T)
    print(li, "rms rho %.4f max|rho| %.3f" % (np.sqrt(np.mean(rho**2)), np.abs(rho).max()),
          " rms T2..T7:", " ".join("%.2e" % np.sqrt(np.mean(T**2)) for T in out),
          " sum T2..T7 rowmean:", " ".join("%.2e" % np.abs(T.sum(1)).mean() for T in out))

    # effect on next-layer pre-activation variance and mean (first order): dvar_i = (W T W^T)_ii
    if li < 15:
        W = Wcol[li+1] if False else None
