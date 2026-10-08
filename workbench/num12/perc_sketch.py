# Percolation (Horvitz-Thompson) sketch of the quenched path-tree term T3 of kappa_3(z_i):
#   T3_i = 3 sum_l W_il a2_l (U_il^2 - V_il),  U = (W*a1) R0,  V_il = sum_k W_ik^2 a1_k^2 R0_kl^2
# Leaves k retained independently with prob p_k, centres l retained with prob q_l; inverse-survival weights.
import numpy as np, sys
from closure import closure, relu_coeffs, phi
n, L, s = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
Ws = list(np.load(f"W_n{n}_L{L}_s{s}.npy")); oc = closure(Ws, keep=True)
rng = np.random.default_rng(1)
for lay in [1, 7, 14]:
    d = oc[lay]; W = Ws[lay+1]; A = relu_coeffs(d["mu"], d["sig"], 3); a1, a2 = A[1], A[2]
    R0 = d["R"].copy(); np.fill_diagonal(R0, 0)
    X = W*a1; U = X @ R0; V = (X*X) @ (R0*R0)
    T3 = 3*np.sum(W*a2*(U*U - V), 1)
    gap = phi(d["mu"]/d["sig"]); neff = gap.sum()**2/(n*(gap**2).sum())
    print(f"layer {lay+1}->{lay+2}: rms T3 {np.sqrt(np.mean(T3**2)):.2e}; effective gap fraction N_eff/n = {neff:.2f}")
    for p in [0.5, 0.25, 0.1]:
        for guide in ["uniform", "gap-guided"]:
            errs = []
            for rep in range(5):
                pk = np.full(n, p)
                if guide == "uniform": ql = np.full(n, p)
                else:
                    w = np.abs(a2); ql = np.minimum(1.0, w/w.mean()*p)    # centre retention proportional to |gap jet|
                Bk = rng.random(n) < pk; Bl = rng.random(n) < ql
                Xs = X[:, Bk]/pk[Bk]; Us = Xs @ R0[Bk][:, Bl]                         # leaves sketched
                corr = ((X[:, Bk]**2)*((1-pk[Bk])/pk[Bk]**2)) @ (R0[Bk][:, Bl]**2)   # e^2=e correction: shared leaf survives once
                Vs = ((X[:, Bk]**2)/pk[Bk]) @ (R0[Bk][:, Bl]**2)
                T3s = 3*np.sum(W[:, Bl]*a2[Bl]/ql[Bl]*(Us*Us - corr - Vs), 1)
                cost = (Bk.mean()*Bl.mean())
                errs.append(np.linalg.norm(T3s-T3)/np.linalg.norm(T3))
            print(f"    p={p:4.2f} {guide:10s}: relative error {np.mean(errs):.2f}  (product cost fraction ~{cost:.3f})")
