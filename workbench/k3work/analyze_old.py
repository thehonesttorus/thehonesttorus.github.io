import numpy as np, pickle, sys
sys.path.insert(0, "../num12")
i = int(sys.argv[1])
mt = np.load(f"../official/truth_off{i}.npz")["m"].astype(np.float64)
pf, pw = np.load(f"pred_off{i}_win0.npy"), np.load(f"pred_off{i}_win4.npy")
df, dw = pickle.load(open(f"dump_off{i}_win0.pkl", "rb")), pickle.load(open(f"dump_off{i}_win4.pkl", "rb"))
# (1) final-mean effect of old sources vs the gain zero mode
from gac import gac
from closure import closure
Ws = list(np.load(f"../official/W_off{i}.npy").astype(np.float64))
mg = gac(Ws)[-1]["m"]; mc = closure(Ws)[-1]["m"]
d_old = pf[-1] - pw[-1]; d_gain = mg - mc; m = pf[-1]
def share(v, u): return (v @ u)**2/((u @ u)*(v @ v))
print(f"final layer: |d_old|^2/n {np.mean(d_old**2):.3e}; share along mean dir {share(d_old, m):.3f}; along GAC gain correction {share(d_old, d_gain):.3f}; corr(d_old, d_gain) {np.corrcoef(d_old, d_gain)[0,1]:+.3f}")
print(f"  errors: full {np.mean((pf[-1]-mt[-1])**2):.3e}, win4 {np.mean((pw[-1]-mt[-1])**2):.3e}, win4 + best scalar*d_gain {min(np.mean((pw[-1]+a*d_gain-mt[-1])**2) for a in np.linspace(-3,3,601)):.3e}, win4 + best scalar*mean {min(np.mean((pw[-1]*(1+a)-mt[-1])**2) for a in np.linspace(-0.01,0.01,401)):.3e}")
# (2) per-layer old-source slices vs the scale-mixture form  kappa3(i,j,k) = g/2 (mu_i S_jk + mu_j S_ik + mu_k S_ij)
print(" layer | |D21_old|/|D21| | R2(D21_old ~ 2 mu_i S_ic + mu_c S_ii) | R2(D21_old ~ a mu_i S_ic + b mu_c S_ii) | |D3_old|/|D3| | R2(D3_old ~ mu_i S_ii) | implied gamma")
for a, b in zip(df, dw):
    if a["D21"] is None: continue
    l = a["layer"]; D21o = a["D21"] - b["D21"]; D3o = a["D3"] - b["D3"]
    if np.abs(D21o).max() == 0: continue
    mu, S = a["mu"], a["C_pre"]; n = len(mu); off = ~np.eye(n, dtype=bool)
    F1 = (mu[:, None]*S)[off]; F2 = (np.diag(S)[None, :]*mu[None, :]*0 + mu[None, :]*np.diag(S)[:, None])[off]
    y = D21o[off]; g = 2*F1 + F2
    al = (g @ y)/(g @ g); r2a = 1 - np.sum((y - al*g)**2)/np.sum(y**2)
    X = np.stack([F1, F2], 1); c, *_ = np.linalg.lstsq(X, y, rcond=None); r2b = 1 - np.sum((y - X @ c)**2)/np.sum(y**2)
    f3 = mu*np.diag(S); be = (f3 @ D3o)/(f3 @ f3); r23 = 1 - np.sum((D3o - be*f3)**2)/np.sum(D3o**2)
    print(f"  {l+1:3d}  |    {np.linalg.norm(D21o)/np.linalg.norm(a['D21'][off]):.3f}       |        {r2a:+.3f}                     |        {r2b:+.3f}                    |   {np.linalg.norm(D3o)/np.linalg.norm(a['D3']):.3f}     |   {r23:+.3f}   | {2*al:+.4f}")
