# Monte Carlo ground truth for He ReLU MLPs (no biases), Gaussian inputs.
import numpy as np, sys, time
n, L, seed, T = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(float(sys.argv[4]))
rng = np.random.default_rng(seed)
Ws = [rng.standard_normal((n, n)) * np.sqrt(2.0/n) for _ in range(L)]
np.save(f"W_n{n}_L{L}_s{seed}.npy", np.array(Ws))
Wf = [W.T.astype(np.float32).copy() for W in Ws]
xr = np.random.default_rng(seed + 10**6)
B = 8192
S1 = np.zeros((L, n)); S2 = np.zeros((L, n)); done = 0
t0 = time.time()
while done < T:
    b = min(B, T - done)
    h = xr.standard_normal((b, n), dtype=np.float32)
    for l in range(L):
        h = np.maximum(h @ Wf[l], 0)
        S1[l] += h.sum(0, dtype=np.float64); S2[l] += (h.astype(np.float64)**2).sum(0)
    done += b
m = S1/T; v = S2/T - m**2
np.savez(f"truth_n{n}_L{L}_s{seed}.npz", m=m, v=v, T=T)
print(n, L, seed, T, "time %.1fs" % (time.time()-t0), "final mean of means %.4f, avg var %.4f" % (m[-1].mean(), v[-1].mean()))
