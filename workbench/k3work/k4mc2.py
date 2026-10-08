# Monte Carlo of the pre-activation (2,2) cumulant slice K22(y_t)_ij = kappa(y_i, y_i, y_j, y_j) at target layers t = s+1,
# together with the post-activation pair cumulants of h_s (C, E[X^2 X'^2], E X^4) at the source layers, two halves.
#   python k4mc2.py NET "3,6,9,12" T
import sys, time, numpy as np
net = int(sys.argv[1]); srcs = [int(v) for v in sys.argv[2].split(",")]; T = int(float(sys.argv[3])); B = 4000
Wc = np.load(f"../official/W_off{net}.npy"); n = Wc.shape[1]
Wt = [np.ascontiguousarray(W.T).astype(np.float32) for W in Wc]
tg = [s + 1 for s in srcs]; top = max(tg); f64 = np.float64
rng = np.random.default_rng(9090 + net); t0 = time.time()
mh = {s: np.zeros(n) for s in srcs}; my = {t: np.zeros(n) for t in tg}; done = 0
while done < 10 ** 6:
    h = rng.standard_normal((B, n), dtype=np.float32)
    for l in range(top + 1):
        y = h @ Wt[l]
        if l in my: my[l] += y.sum(0, dtype=f64)
        h = np.maximum(y, 0)
        if l in mh: mh[l] += h.sum(0, dtype=f64)
    done += B
mh = {k: (v / done).astype(np.float32) for k, v in mh.items()}; my = {k: (v / done).astype(np.float32) for k, v in my.items()}
def blank(): return dict(C=np.zeros((n, n)), M22=np.zeros((n, n)), s1=np.zeros(n), s2=np.zeros(n), s4=np.zeros(n))
acc = [{**{("h", s): blank() for s in srcs}, **{("y", t): blank() for t in tg}} for _ in range(2)]
def add(a, X):
    X2 = X * X
    a["C"] += (X.T @ X).astype(f64); a["M22"] += (X2.T @ X2).astype(f64)
    a["s1"] += X.sum(0, dtype=f64); a["s2"] += X2.sum(0, dtype=f64); a["s4"] += (X2 * X2).sum(0, dtype=f64)
done = 0; k = 0
while done < T:
    a2 = acc[k % 2]; h = rng.standard_normal((B, n), dtype=np.float32)
    for l in range(top + 1):
        y = h @ Wt[l]
        if l in my: add(a2[("y", l)], y - my[l])
        h = np.maximum(y, 0)
        if l in mh: add(a2[("h", l)], h - mh[l])
    done += B; k += 1
    if k % 100 == 0: print(f"  {done}/{T} {time.time()-t0:.0f}s", flush=True)
np.savez(f"k4mc2_off{net}.npz", T=done // 2, **{f"h{hf}_{kind}{l}_{key}": v for hf in range(2) for (kind, l), d in acc[hf].items() for key, v in d.items()})
print(f"done {done} {time.time()-t0:.0f}s", flush=True)
