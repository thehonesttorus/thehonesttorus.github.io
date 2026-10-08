# Monte Carlo first and second moments of the post-activations h_l at every layer (for the local-defect ledger).
import numpy as np, sys
n, L, s, T = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(float(sys.argv[4]))
Ws = np.load(f"W_n{n}_L{L}_s{s}.npy"); Wf = [W.T.astype(np.float32).copy() for W in Ws]
rng = np.random.default_rng(s + 31337); B = 8192
S1 = np.zeros((L, n)); S2 = np.zeros((L, n, n)); done = 0
while done < T:
    b = min(B, T-done); h = rng.standard_normal((b, n), dtype=np.float32)
    for l in range(L):
        h = np.maximum(h @ Wf[l], 0)
        S1[l] += h.sum(0, dtype=np.float64); S2[l] += (h.T @ h).astype(np.float64)
    done += b
np.savez(f"hmom_n{n}_L{L}_s{s}.npz", S1=S1/T, S2=S2/T, T=T)
print("done", T)
