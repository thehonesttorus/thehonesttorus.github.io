# Old-source (2,1) slice vs the "gain field" form  kappa3(i,j,k) = v_i S_jk + v_j S_ik + v_k S_ij  (v in R^n, fitted)
import numpy as np, pickle, sys
i = int(sys.argv[1]); df, dw = pickle.load(open(f"dump_off{i}_win0.pkl", "rb")), pickle.load(open(f"dump_off{i}_win4.pkl", "rb"))
print(" layer | R2 scalar-gain form | R2 gain-field form (v fitted, n params) | R2 field on ALL of D21 (not just old) | corr(v, mu)")
for a, b in zip(df, dw):
    if a["D21"] is None: continue
    l = a["layer"]; D21o = a["D21"] - b["D21"]
    if np.abs(D21o).max() == 0: continue
    mu, S = a["mu"], a["C_pre"]; n = len(mu); d = np.diag(S).copy(); So = S.copy(); np.fill_diagonal(So, 0)
    def fit(D):
        D = D.copy(); np.fill_diagonal(D, 0)
        A = 2*So*(d[:, None] + d[None, :]); np.fill_diagonal(A, 4*(So**2).sum(1) + (d**2).sum() - d**2)
        bvec = 2*(So*D).sum(1) + (d[:, None]*D).sum(0)
        v = np.linalg.solve(A, bvec)
        P = 2*v[:, None]*So + d[:, None]*v[None, :]; np.fill_diagonal(P, 0)
        return v, 1 - np.sum((D - P)**2)/np.sum(D**2)
    g = 2*mu[:, None]*So + mu[None, :]*d[:, None]; np.fill_diagonal(g, 0); Do = D21o.copy(); np.fill_diagonal(Do, 0)
    al = np.sum(g*Do)/np.sum(g*g); r2s = 1 - np.sum((Do - al*g)**2)/np.sum(Do**2)
    v, r2v = fit(D21o); _, r2all = fit(a["D21"])
    print(f"  {l+1:3d}  |      {r2s:+.3f}         |              {r2v:+.3f}                    |            {r2all:+.3f}                 | {np.corrcoef(v, mu)[0,1]:+.3f}")
