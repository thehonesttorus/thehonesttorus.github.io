# Monte Carlo of the pair cumulant slices needed for the two-loop gain derivation, official network NET.
# Per layer l: pre-activation y_l (= z) and post-activation h_l = relu(y_l):
#   y: E y, E y y', E y^2 y' (n x n), E y^3 y' (n x n), E y^2 y'^2 (n x n), E y^3, E y^4   -> kappa3 slices D3, D21 and kappa4 slices
#   h: E h, E h h', E h^2 h', E h^2 h'^2, E h^3, E h^4                                      -> true pair cumulants C3, C4 of h
#   norms: |h|^2 and |z|^2 first two moments (measured gamma per layer)
import numpy as np, sys, time
net = int(sys.argv[1]) if len(sys.argv) > 1 else 0
T = int(float(sys.argv[2])) if len(sys.argv) > 2 else 400_000
Wc = np.load(f"../official/W_off{net}.npy"); n, L = Wc.shape[1], Wc.shape[0]
Wf = [np.ascontiguousarray(W.T).astype(np.float32) for W in Wc]
rng = np.random.default_rng(777)
f64 = np.float64
acc = {k: np.zeros((L, n, n), f64) for k in ("S2y", "M21y", "M31y", "M22y", "S2h", "M21h", "M22h")}
vec = {k: np.zeros((L, n), f64) for k in ("s1y", "s3y", "s4y", "s1h", "s3h", "s4h")}
nrm = np.zeros((L, 4), f64)   # qs1, qs2 (|h|^2), rs1, rs2 (|z|^2)
B = 4096; done = 0; t0 = time.time()
while done < T:
    b = min(B, T - done)
    h = rng.standard_normal((b, n), dtype=np.float32)
    for l in range(L):
        z = h @ Wf[l]; z2 = z * z; z3 = z2 * z
        vec["s1y"][l] += z.sum(0); vec["s3y"][l] += z3.sum(0); vec["s4y"][l] += (z2 * z2).sum(0)
        acc["S2y"][l] += z.T @ z; acc["M21y"][l] += z2.T @ z; acc["M31y"][l] += z3.T @ z; acc["M22y"][l] += z2.T @ z2
        h = np.maximum(z, 0); h2 = h * h
        vec["s1h"][l] += h.sum(0); vec["s3h"][l] += (h2 * h).sum(0); vec["s4h"][l] += (h2 * h2).sum(0)
        acc["S2h"][l] += h.T @ h; acc["M21h"][l] += h2.T @ h; acc["M22h"][l] += h2.T @ h2
        q = h2.sum(1).astype(f64); r = z2.sum(1).astype(f64)
        nrm[l] += [q.sum(), (q * q).sum(), r.sum(), (r * r).sum()]
    done += b
    if done % (B * 10) == 0: print(f"{done} samples, {time.time()-t0:.0f}s", flush=True)
out = {k: v / T for k, v in acc.items()}; out.update({k: v / T for k, v in vec.items()}); out["nrm"] = nrm / T; out["T"] = T
np.savez(f"mc_cum_off{net}.npz", **out)
print(f"done: {T} samples in {time.time()-t0:.0f}s", flush=True)
