# Is the gate-covariance correction the shadow of a small separator?  Recompute D3corr(l+1) on the same output neurons as
# gatecov.py with the pre-activation correlation rho replaced by its rank-K truncation (top-K eigen-directions of the
# off-diagonal correlation, diagonal removed), K = 1, 4, 16, 64, 256, and report the fraction of the exact correction
# captured (correlation, slope, residual rms).  A small hidden variable screening the gate correlations would show up as a
# large captured fraction at small K.
import numpy as np, pickle, sys, time
from scipy.special import ndtr
net = int(sys.argv[1]) if len(sys.argv) > 1 else 0; layer = int(sys.argv[2]) if len(sys.argv) > 2 else 10; NI = int(sys.argv[3]) if len(sys.argv) > 3 else 48
KS = [1, 4, 16, 64, 256]
P = pickle.load(open(f"v29legs_off{net}.pkl", "rb")); legs = {d["layer"]: d for d in P["legs"]}; dumps = {d["layer"]: d for d in P["dumps"]}
Lg = legs[layer]; d = dumps[layer]
mu, var, C_off = d["mu"], d["var"], d["C_off"]; n = len(mu); sig = np.sqrt(var); a = mu / sig; Ph = ndtr(a); ph = np.exp(-a*a/2)/np.sqrt(2*np.pi)
R = (C_off / np.outer(sig, sig)).astype(np.float64)
Wcol = np.load(f"../official/W_off{net}.npy"); Wn = Wcol[layer + 1].astype(np.float64)
A, Pl, Z, L = Lg["A"].astype(np.float64), Lg["P"].astype(np.float64), Lg["Z"].astype(np.float64), Lg["L"].astype(np.float64)
k = A.shape[0]; w2b = [np.asarray(x, dtype=np.float64) for x in Lg["w2b"]]; s_l = [np.asarray(x, dtype=np.float64) for x in Lg["s"]]; e_l = [np.asarray(x, dtype=np.float64) for x in Lg["e"]]
M = [Pl[s] * s_l[s][None, :] + 3.0 * A[s] * e_l[s][None, :] + Z[s] @ L[s].T for s in range(k)]
if "w1" in d and "S3c" in d:
    A_new = d["w1"][:, None] * C_off; P_new = np.eye(n); M_new = np.diag(d["S3c"]) + 3.0 * A_new * d["e_b"][None, :] + 3.0 * d["Rres"]
    A = np.concatenate([A, A_new[None]], 0); Pl = np.concatenate([Pl, P_new[None]], 0); M.append(M_new); w2b.append(np.asarray(d["w2"], dtype=np.float64)); k += 1
# eigen-directions of the off-diagonal correlation (symmetric, zero diagonal)
ev, U = np.linalg.eigh(R); order = np.argsort(-np.abs(ev)); ev, U = ev[order], U[:, order]
print(f"net {net} layer {layer}: {k} sources; |eigenvalues| of the off-diagonal correlation: top {np.abs(ev[:4]).round(2)}, 16th {abs(ev[15]):.2f}, 64th {abs(ev[63]):.2f}, 256th {abs(ev[255]):.2f}; sum of squares captured by top K: " + ", ".join(f"K={K}: {np.sum(ev[:K]**2)/np.sum(ev**2):.3f}" for K in KS))
rng = np.random.default_rng(0); idx = np.sort(rng.choice(n, NI, replace=False))
def contraction(Rm):
    out = np.zeros(NI)
    for t, i in enumerate(idx):
        w = Wn[i]; v = w * Ph; u = w * ph; K = (u[:, None] * Rm) * u[None, :]; cr = 0.0
        for s in range(k):
            As, Ps, Ms = A[s], Pl[s], M[s]; Av, Pv, Mv = As.T @ v, Ps.T @ v, Ms.T @ v
            KA = K @ As; KP = K @ Ps
            dAA = np.sum(As * KA, axis=0); dAP = np.sum(Ps * KA, axis=0); dPP = np.sum(Ps * KP, axis=0); dMP = np.sum(Ms * KP, axis=0)
            cr += 3.0 * (np.sum(w2b[s] * (dAA * Pv + 2.0 * dAP * Av)) + (1.0 / 3.0) * np.sum(2.0 * dMP * Pv + dPP * Mv))
        out[t] = cr
    return out
t0 = time.time(); exact = contraction(R); rms = lambda x: np.sqrt(np.mean(x**2)); print(f"exact correction: rms {rms(exact):.3e}  [{time.time()-t0:.0f}s]", flush=True)
for K in KS:
    RK = (U[:, :K] * ev[:K]) @ U[:, :K].T; np.fill_diagonal(RK, 0.0)
    cK = contraction(RK); res = exact - cK
    print(f"  rank {K:3d}: captured corr {np.corrcoef(cK, exact)[0,1]:+.3f}, slope {np.sum(cK*exact)/np.sum(exact**2):.3f}, residual rms {rms(res):.3e} = {rms(res)/rms(exact):.1%} of exact  [{time.time()-t0:.0f}s]", flush=True)
np.savez(f"gatecov_rank_off{net}_l{layer}.npz", idx=idx, exact=exact, ev=ev)
