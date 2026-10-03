# Is the residue the gain (scale-mixture) variance?  Compare the oracle scale coefficient c_l of the closure
# error with the coherent normalized fourth cumulant gamma_l = avg_{i!=j} k(z_i,z_i,z_j,z_j)/(s_i^2 s_j^2).
import numpy as np, sys
from closure import closure
n, L, s, T = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(float(sys.argv[4]))
Ws = list(np.load(f"W_n{n}_L{L}_s{s}.npy")); tr = np.load(f"truth_n{n}_L{L}_s{s}.npz")
oc = closure(Ws); Wf = [W.T.astype(np.float32) for W in Ws]
rng = np.random.default_rng(5); B = 10000
S1 = np.zeros((L, n)); S2 = np.zeros((L, n)); Q = np.zeros((L, n, n)); G2 = np.zeros((L, n, n)); N2 = np.zeros(L); N4 = np.zeros(L)
for _ in range(T//B):
    h = rng.standard_normal((B, n)).astype(np.float32)
    for l in range(L):
        z = (h @ Wf[l]).astype(np.float64); h = np.maximum(z, 0).astype(np.float32)
        S1[l] += z.sum(0); S2[l] += (z*z).sum(0); Q[l] += z.T @ z
        G2[l] += (z*z).T @ (z*z)
        nn = (h.astype(np.float64)**2).sum(1); N2[l] += nn.sum(); N4[l] += (nn*nn).sum()
for l in range(L):
    mu = S1[l]/T; M2 = Q[l]/T; C = M2 - np.outer(mu, mu)
    # raw moment E[z_i^2 z_j^2]; for near-zero-mean approximation use central-ish version via Gaussian subtraction
    E22 = G2[l]/T; Ez2 = S2[l]/T
    k22 = E22 - np.outer(Ez2, Ez2) - 2*M2**2 + 2*np.outer(mu**2, mu**2)   # Gaussian (Isserlis) part removed
    off = ~np.eye(n, dtype=bool)
    gam = np.mean((k22/np.outer(Ez2, Ez2))[off])
    m = tr["m"][l]; mc = oc[l]["m"]; c = mc @ (mc - m)/(mc @ mc)
    gvar = N4[l]/T/(N2[l]/T)**2 - 1
    print(f"layer {l+1:2d}: oracle scale c {c:+.4f} | gamma (coherent K22 / s^2 s^2) {gam:+.4f} -> gamma/8 {gam/8:+.4f} | Var(|h|^2)/E^2 {gvar:.4f}")
