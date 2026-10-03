# Degree-resolved variance shares on the sphere, measured with designs (no truth needed):
#   one random orthonormal frame (n pts) kills degree 2 exactly  -> Var_frame/Var_iid = 1 - v2 (other degrees ~iid)
#   antipodal pair average keeps only even degrees              -> Var_pair = even share (per pair, vs Var/2)
import numpy as np, sys
def net(n, L, seed):
    rng = np.random.default_rng(seed); return [(rng.standard_normal((n, n))*np.sqrt(2/n)).astype(np.float32) for _ in range(L)]
def fwd(Wf, X):
    h = X
    for W in Wf: h = np.maximum(h @ W, 0)
    return h.astype(np.float64)
def shares(Wf, n, reps=96, seed=0):
    rng = np.random.default_rng(seed); fr, iid, pr = [], [], []
    for _ in range(reps):
        Q = np.linalg.qr(rng.standard_normal((n, n)))[0].astype(np.float32)
        fr.append(fwd(Wf, Q).mean(0))                                  # one frame on the unit sphere
        T = rng.standard_normal((n, n)).astype(np.float32); T /= np.linalg.norm(T, axis=1, keepdims=True)
        iid.append(fwd(Wf, T).mean(0))                                 # n iid points on the unit sphere
        U = T[: n//2]; pr.append((fwd(Wf, U).mean(0) + fwd(Wf, -U).mean(0))/2)   # n/2 antipodal pairs
    vf, vi, vp = np.var(fr, 0).mean(), np.var(iid, 0).mean(), np.var(pr, 0).mean()
    return 1 - vf/vi, vp/vi   # degree-2 share; pair variance per pass relative to iid (= 2 x even share)
for (n, L, src) in [(1024, 16, "official"), (256, 32, "he"), (256, 16, "he"), (1024, 32, "he")]:
    if src == "official":
        Wc = np.load("../official/W_off0.npy"); Wf = [np.ascontiguousarray(W.T) for W in Wc]
    else: Wf = net(n, L, 7)
    v2, pair = shares(Wf, n, reps=64 if n == 1024 else 160)
    print(f"n={n:5d} L={L:2d}: degree-2 share {v2:.3f}; even share {pair/2:.3f} (odd {1-pair/2:.3f}); predicted frame+antipode factor 2*(even - v2) = {2*(pair/2 - v2):.3f}", flush=True)
