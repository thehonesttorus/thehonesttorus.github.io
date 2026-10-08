# Monte Carlo for the per-neuron fourth cumulant of z_(s+1) = W_(s+1) h_s by index class (note XXI's power sums), plus the
# true post-activation pair cumulants of h_s (C, E[X^2 X'^2], E X^4) so the pair class can be rebuilt by any contraction.
# Two halves (even / odd batches) for a noise estimate.   python k4mc.py NET "3,6,9,12" T
import sys, time, numpy as np
net = int(sys.argv[1]); srcs = [int(v) for v in sys.argv[2].split(",")]; T = int(float(sys.argv[3])); B = 4000
Wc = np.load(f"../official/W_off{net}.npy"); L, n = Wc.shape[0], Wc.shape[1]
Wt = [np.ascontiguousarray(W.T).astype(np.float32) for W in Wc]
Wp = {s: [np.ascontiguousarray((Wc[s + 1].astype(np.float64) ** r).T).astype(np.float32) for r in (1, 2, 3, 4)] for s in srcs}
top = max(srcs) + 1; f64 = np.float64
rng = np.random.default_rng(4242 + net); t0 = time.time()
m = {s: np.zeros(n) for s in srcs}; done = 0
while done < 10 ** 6:                                    # pre-pass: post-activation means for centring
    h = rng.standard_normal((B, n), dtype=np.float32)
    for l in range(top):
        h = np.maximum(h @ Wt[l], 0)
        if l in srcs: m[l] += h.sum(0, dtype=f64)
    done += B
m = {s: (m[s] / done).astype(np.float32) for s in srcs}
print(f"pre-pass {done} samples {time.time()-t0:.0f}s", flush=True)
acc = [{s: dict(C=np.zeros((n, n)), M22=np.zeros((n, n)), s1=np.zeros(n), s2=np.zeros(n), s4=np.zeros(n),
                E4=np.zeros(n), E13=np.zeros(n), E22=np.zeros(n), E112=np.zeros(n), E1111=np.zeros(n), E11=np.zeros(n),
                z1=np.zeros(n), z2=np.zeros(n), z3=np.zeros(n), z4=np.zeros(n)) for s in srcs} for _ in range(2)]
done = 0; k = 0
while done < T:
    a2 = acc[k % 2]
    h = rng.standard_normal((B, n), dtype=np.float32)
    for l in range(top + 1):
        y = h @ Wt[l]
        if l - 1 in srcs:                                # raw moments of the target pre-activation z_(s+1)
            a = a2[l - 1]; y64 = y.astype(f64); y2 = y64 * y64
            a["z1"] += y64.sum(0); a["z2"] += y2.sum(0); a["z3"] += (y2 * y64).sum(0); a["z4"] += (y2 * y2).sum(0)
        if l == top: break
        h = np.maximum(y, 0)
        if l in srcs:
            a = a2[l]; X = h - m[l]; X2 = X * X
            a["C"] += (X.T @ X).astype(f64); a["M22"] += (X2.T @ X2).astype(f64)
            a["s1"] += X.sum(0, dtype=f64); a["s2"] += X2.sum(0, dtype=f64); a["s4"] += (X2 * X2).sum(0, dtype=f64)
            P1 = (X @ Wp[l][0]).astype(f64); P2 = (X2 @ Wp[l][1]).astype(f64)
            P3 = ((X2 * X) @ Wp[l][2]).astype(f64); P4 = ((X2 * X2) @ Wp[l][3]).astype(f64)
            a["E4"] += P4.sum(0); a["E13"] += (P1 * P3).sum(0); a["E22"] += (P2 * P2).sum(0)
            a["E112"] += (P1 * P1 * P2).sum(0); a["E1111"] += (P1 ** 4).sum(0); a["E11"] += (P1 * P1).sum(0)
    done += B; k += 1
    if k % 100 == 0: print(f"  {done}/{T} {time.time()-t0:.0f}s", flush=True)
np.savez(f"k4mc_off{net}.npz", T=done // 2, srcs=np.array(srcs),
         **{f"h{hf}_s{s}_{key}": v for hf in range(2) for s in srcs for key, v in acc[hf][s].items()})
print(f"done {done} samples {time.time()-t0:.0f}s", flush=True)
