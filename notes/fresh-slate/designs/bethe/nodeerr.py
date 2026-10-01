import sys, numpy as np
sys.path.insert(0, '../../bench'); sys.path.insert(0, '.')
import bench, bethe, oracle
S = bench.load_set(sys.argv[1]); i = int(sys.argv[2]); N = float(sys.argv[3])
W = bench.weights(S, i).astype(np.float64)
tr = oracle.mc_state(W, int(N))
# run own chain, record node states
Ls, n, _ = W.shape
m = np.zeros(n); C = W[0].T @ W[0]; K = np.zeros((n, n)); k4 = np.zeros(n); prev = None
rel = lambda a, b: np.sqrt(np.mean((a - b) ** 2)) / np.sqrt(np.mean(b ** 2))
print('layer | rel err m, v, k3, k4 | rms true k3, k4 | corr(k4 own,true)')
for l in range(Ls):
    t = tr[l]
    print(f"{l:2d} | {rel(m,t['m']):.3f} {rel(np.diag(C),np.diag(t['C'])):.4f} {rel(np.diag(K),np.diag(t['K'])):.3f} {rel(k4,t['k4']):.3f} | "
          f"{np.sqrt(np.mean(np.diag(t['K'])**2)):.4f} {np.sqrt(np.mean(t['k4']**2)):.4f} | {np.corrcoef(k4,t['k4'])[0,1] if l else 0:.3f} mean k4 own {k4.mean():.4f} true {t['k4'].mean():.4f}", flush=True)
    mu, Ca, Ka, k3a, k4a, c, L0, Lm1 = bethe.relu_map_edges(m, C, K, k4)
    if l + 1 == Ls: break
    Wn = W[l + 1]
    M = c * c * (L0 ** 2)[:, None] * Lm1[None, :]
    spec = (k3a, Ka - M, c, L0, Lm1)
    Kn = bethe.contract(*spec, Wn, Wn)
    if prev is not None:
        pspec, Wl = prev
        P = Wl @ (L0[:, None] * Wn)
        Kz = K.copy(); k3z = np.diag(K).copy(); np.fill_diagonal(Kz, 0.0)
        Kn = Kn + bethe.contract(*pspec, P, P) - bethe.contract(k3z * L0 ** 3, Kz * (L0 ** 2)[:, None] * L0[None, :], None, L0, Lm1, Wn, Wn)
    prev = (spec, Wn)
    m = mu @ Wn; C = Wn.T @ Ca @ Wn; K = Kn; k4 = (Wn ** 4).T @ k4a
