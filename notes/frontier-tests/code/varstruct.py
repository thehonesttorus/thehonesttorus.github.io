# Variance structure of final-layer activations on an official network: antithetic, common mode, top PCs, previous layer.
import numpy as np, time
Wcol = np.load("W_off0.npy"); Wf = [np.ascontiguousarray(W.T) for W in Wcol]   # h @ Wf = relu(W h) row form
rng = np.random.default_rng(0); N = 16384; B = 4096
H, Ha, Hp = [], [], []
t0 = time.time()
for b in range(N//B):
    x = rng.standard_normal((B, 1024)).astype(np.float32)
    for sgn, store in [(1, H), (-1, Ha)]:
        h = sgn*x
        for l, W in enumerate(Wf):
            if l == 15 and sgn == 1: Hp.append(h.copy())
            h = np.maximum(h @ W, 0)
        store.append(h)
H = np.concatenate(H).astype(np.float64); Ha = np.concatenate(Ha).astype(np.float64); Hp = np.concatenate(Hp).astype(np.float64)
print(f"forward {time.time()-t0:.0f}s")
v = H.var(0); print(f"avg per-neuron var {v.mean():.4f}, avg mean {H.mean():.4f}")
va = ((H + Ha)/2).var(0); print(f"antithetic pair-average var / (var/2) = {va.mean()/(v.mean()/2):.3f}  (1 = no gain)")
s = H.mean(1); sc = s - s.mean(); beta = (H - H.mean(0)).T @ sc / (sc @ sc)
r = (H - H.mean(0)) - np.outer(sc, beta); print(f"common mode s(x)=mean_j h_j: explains {1 - r.var(0).mean()/v.mean():.3f} of variance")
C = np.cov(H.T); ev, V = np.linalg.eigh(C); ev = ev[::-1]
print("top PC variance shares:", " ".join(f"{e/ev.sum():.3f}" for e in ev[:8]), f"| top-1..32 cumulative {ev[:32].sum()/ev.sum():.3f}, top-128 {ev[:128].sum()/ev.sum():.3f}")
# linear predictor from previous layer (exact relation h_L = relu(W h_{L-1}))
X = Hp - Hp.mean(0); Y = H - H.mean(0)
coef, *_ = np.linalg.lstsq(X, Y, rcond=None); res = Y - X @ coef
print(f"linear in h_15 explains {1 - res.var(0).mean()/v.mean():.3f} (in-sample, p/N={1024/N:.3f})")
