# The gate-covariance correction to the all-distinct third-cumulant transport, computed exactly for one layer from the
# chain's live source legs (dumped by run_v29legs.py), on a subset of output neurons of the next layer.
#   T_abc(l) = sum_j w2_j Sym(A_aj A_bj P_cj) + (1/3) sum_j Sym(M_aj P_bj P_cj),  M = P d(s) + 3 A d(e) + Z L^T  (per source)
#   product gates (chain):  D3prod_i(l+1) = sum_abc W_ia W_ib W_ic Phi_a Phi_b Phi_c T_abc
#   first-order correction: D3corr_i(l+1) = 3 sum_abc W_ia W_ib W_ic rho_ab phi_a phi_b Phi_c T_abc   (rho = pre-activation correlation, a != b)
import numpy as np, pickle, sys, time
from scipy.special import ndtr
net = int(sys.argv[1]) if len(sys.argv) > 1 else 0; layer = int(sys.argv[2]) if len(sys.argv) > 2 else 10; NI = int(sys.argv[3]) if len(sys.argv) > 3 else 48
P = pickle.load(open(f"v29legs_off{net}.pkl", "rb")); legs = {d["layer"]: d for d in P["legs"]}; dumps = {d["layer"]: d for d in P["dumps"]}
Lg = legs[layer]; d = dumps[layer]
mu, var, C_off = d["mu"], d["var"], d["C_off"]; n = len(mu); sig = np.sqrt(var); a = mu / sig; Ph = ndtr(a); ph = np.exp(-a*a/2)/np.sqrt(2*np.pi)
R = C_off / np.outer(sig, sig)   # zero diagonal: the correction is over a != b
Wcol = np.load(f"../official/W_off{net}.npy"); Wn = Wcol[layer + 1].astype(np.float64)   # the chain's W at layer+1 (W = w32.T, w = stored.T)
A, Pl, Z, L = Lg["A"].astype(np.float64), Lg["P"].astype(np.float64), Lg["Z"].astype(np.float64), Lg["L"].astype(np.float64)
k = A.shape[0]; w2b = [np.asarray(x, dtype=np.float64) for x in Lg["w2b"]]; s_l = [np.asarray(x, dtype=np.float64) for x in Lg["s"]]; e_l = [np.asarray(x, dtype=np.float64) for x in Lg["e"]]
print(f"net {net} layer {layer}: {k} sources, legs {A.shape}, Z {Z.shape}; chain D3(l) rms {np.sqrt(np.mean(Lg['D3']**2)):.3e}")
M = [Pl[s] * s_l[s][None, :] + 3.0 * A[s] * e_l[s][None, :] + Z[s] @ L[s].T for s in range(k)]
# the source born at this layer (not yet in the stacks at the contraction): A = d(w1) C_off, P = I, M = d(S3c) + 3 A d(e_b) + 3 Rres
if "w1" in d and "S3c" in d:
    A_new = d["w1"][:, None] * C_off; P_new = np.eye(n); M_new = np.diag(d["S3c"]) + 3.0 * A_new * d["e_b"][None, :] + 3.0 * d["Rres"]
    A = np.concatenate([A, A_new[None]], 0); Pl = np.concatenate([Pl, P_new[None]], 0); M.append(M_new)
    w2b.append(np.asarray(d["w2"], dtype=np.float64)); k += 1; print(f"  newborn of layer {layer} added as source {k}")
rng = np.random.default_rng(0); idx = np.sort(rng.choice(n, NI, replace=False))
prod = np.zeros(NI); corr = np.zeros(NI); t0 = time.time()
for t, i in enumerate(idx):
    w = Wn[i]; v = w * Ph; u = w * ph
    K = (u[:, None] * R) * u[None, :]            # K_ab = W_ia W_ib rho_ab phi_a phi_b, zero diagonal
    pr = 0.0; cr = 0.0
    for s in range(k):
        As, Ps, Ms = A[s], Pl[s], M[s]
        Av, Pv, Mv = As.T @ v, Ps.T @ v, Ms.T @ v
        pr += 3.0 * np.sum(w2b[s] * Av * Av * Pv) + np.sum(Mv * Pv * Pv)
        KA = K @ As; KP = K @ Ps
        dAA = np.sum(As * KA, axis=0); dAP = np.sum(Ps * KA, axis=0); dPP = np.sum(Ps * KP, axis=0); dMP = np.sum(Ms * KP, axis=0)
        cr += 3.0 * (np.sum(w2b[s] * (dAA * Pv + 2.0 * dAP * Av)) + (1.0 / 3.0) * np.sum(2.0 * dMP * Pv + dPP * Mv))
    prod[t] = pr; corr[t] = cr
    if t % 8 == 7: print(f"  {t+1}/{NI} neurons, {time.time()-t0:.0f}s", flush=True)
D3next = dumps[layer + 1]["D3"][idx] if (layer + 1) in dumps else None
rms = lambda x: np.sqrt(np.mean(x**2))
print(f"product-gate D3(l+1) from the legs: rms {rms(prod):.3e}" + (f"; chain's own D3(l+1) on these neurons: rms {rms(D3next):.3e}, corr {np.corrcoef(prod, D3next)[0,1]:+.3f}, slope {np.sum(prod*D3next)/np.sum(D3next**2):.3f}" if D3next is not None else ""))
print(f"gate-covariance correction D3corr(l+1): rms {rms(corr):.3e} = {rms(corr)/rms(prod):.2%} of the product-gate value; corr(D3corr, D3prod) {np.corrcoef(corr, prod)[0,1]:+.3f}; mean/rms {corr.mean()/rms(corr):+.2f}")
# implied mean error at layer l+1: delta m = D3corr/6 * E f''' with E f''' = -alpha phi / sigma^2 at layer l+1
d1 = dumps[layer + 1]; mu1, var1 = d1["mu"][idx], d1["var"][idx]; s1 = np.sqrt(var1); a1 = mu1 / s1; ph1 = np.exp(-a1*a1/2)/np.sqrt(2*np.pi)
dm = corr / 6.0 * (-a1 * ph1 / var1)
print(f"implied post-activation mean error at layer {layer+1} from the dropped correction: rms {rms(dm):.3e} (the output error is 1.48e-4 rms; a per-layer injection of 3-4e-5 rms accumulating in quadrature over the depth reproduces it)")
np.savez(f"gatecov_off{net}_l{layer}.npz", idx=idx, prod=prod, corr=corr, dm=dm)
