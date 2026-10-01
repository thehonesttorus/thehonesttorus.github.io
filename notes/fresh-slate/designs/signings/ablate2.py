import sys, numpy as np, copula1
w, seed = int(sys.argv[1]), int(sys.argv[2])
W = np.load(f"/root/sg/truth/W_w{w}_s{seed}.npy").astype(float)
O = dict(np.load(f"/root/sg/layerstats_w{w}_s{seed}.npz"))
r = lambda x: np.sqrt(np.mean(x ** 2))
n = W.shape[1]; off = ~np.eye(n, dtype=bool)
for use in [('marg', 'cov', 'D')]:
    est, dg = copula1.estimate(W, {'oracle': dict(O, use=use)}, return_diag=True)
    for l in [2, 4, 8, 12, 14]:   # dg[l] describes z_{l+2} (index l+1) from forced state at index l
        d = dg[l]; t = l + 1
        corr = O['Cz'][t] / np.sqrt(np.outer(O['var'][t], O['var'][t]))
        print(f"z idx {t}: var rel {r(d['vz']/O['var'][t]-1):.1e}  k3 {r(d['k3']-O['k3'][t]):.1e}/{r(O['k3'][t]):.1e}  k4 {r(d['k4']-O['k4'][t]):.1e}/{r(O['k4'][t]):.1e}  offcov {r((d['Cz']-O['Cz'][t])[off]):.1e}/{r(O['Cz'][t][off]):.1e}  |corr| mean {np.abs(corr[off]).mean():.2f}  mean corr {corr[off].mean():.2f}")
