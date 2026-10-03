# MC statistics of pre-activations z (mean, var, k3, k4) per layer
import numpy as np, sys
n, L, s, T = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(float(sys.argv[4]))
Ws = np.load(f"W_n{n}_L{L}_s{s}.npy"); Wf = [W.T.astype(np.float32).copy() for W in Ws]
xr = np.random.default_rng(s + 7*10**6); B = 8192
M = np.zeros((4, L, n)); done = 0
while done < T:
    b = min(B, T-done); h = xr.standard_normal((b, n), dtype=np.float32)
    for l in range(L):
        z = h @ Wf[l]; zd = z.astype(np.float64)
        for p in range(4): M[p, l] += (zd**(p+1)).sum(0)
        h = np.maximum(z, 0)
    done += b
M /= T
mu = M[0]; c2 = M[1]-mu**2; c3 = M[2]-3*mu*M[1]+2*mu**3
c4 = M[3]-4*mu*M[2]+6*mu**2*M[1]-3*mu**4 - 3*c2**2
np.savez(f"zstats_n{n}_L{L}_s{s}.npz", mu=mu, var=c2, k3=c3, k4=c4, T=T)
print("done")
