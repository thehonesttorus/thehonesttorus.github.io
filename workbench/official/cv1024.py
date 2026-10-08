# Width-1024 check of the "free" layer-1 control variates and the antithetic pair (official network 0).
import numpy as np
Wcol = np.load("W_off0.npy").astype(np.float64); Wf = [np.ascontiguousarray(W.T).astype(np.float32) for W in Wcol]
mu1 = np.sqrt((Wcol[0]**2).sum(1))/np.sqrt(2*np.pi)                          # exact layer-1 means
U, s, Vt = np.linalg.svd(Wcol[0]); rng = np.random.default_rng(1); N = 32768; B = 4096
H1, HL, X = [], [], []
for b in range(N//B):
    x = rng.standard_normal((B, 1024)).astype(np.float32); h = x
    for l, W in enumerate(Wf):
        h = np.maximum(h @ W, 0)
        if l == 0: H1.append(h.astype(np.float64))
    HL.append(h.astype(np.float64)); X.append(x.astype(np.float64))
H1 = np.concatenate(H1) - mu1; HL = np.concatenate(HL); X = np.concatenate(X); Y = HL - HL.mean(0); v = Y.var(0).mean()
def r2(F):
    F = F - F.mean(0); coef, *_ = np.linalg.lstsq(F, Y, rcond=None); return 1 - (Y - F @ coef).var(0).mean()/v
for k in [8, 64, 256]:
    print(f"top-{k} left-singular coords of h1-mu: R^2 {r2(H1 @ U[:, :k]):.3f}  (in-sample, p/N={k/N:.4f})")
F = np.concatenate([X @ Vt[:2].T, (X @ Vt[:2].T)**2 - 1, H1 @ U[:, :8]], 1); print(f"theory directions + Hermite-2 + top-8 h1 coords: R^2 {r2(F):.3f}")
