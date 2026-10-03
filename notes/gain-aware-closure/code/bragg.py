# Test of the "Bragg peaks" claim: does the law of the deep representation have a point-spectrum component?
# Empirical characteristic function |E exp(i t v.z)| of final-layer pre-activations along random directions,
# along single neurons, and of the final code bits' sign field (Walsh concentration of the code law).
import numpy as np
n, L, s, T = 256, 16, 0, 200000
Ws = np.load(f"W_n{n}_L{L}_s{s}.npy"); Wf = [W.T.astype(np.float32) for W in Ws]
rng = np.random.default_rng(9); Z = []
for _ in range(T//10000):
    h = rng.standard_normal((10000, n)).astype(np.float32)
    for W in Wf[:-1]: h = np.maximum(h @ W, 0)
    Z.append((h @ Wf[-1]).astype(np.float64))
Z = np.concatenate(Z)
ts = np.array([0.5, 1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64])
print(f"noise floor ~ {1/np.sqrt(T):.4f}")
for name, proj in [("random direction 1", rng.standard_normal(n)), ("random direction 2", rng.standard_normal(n)),
                   ("top covariance eigvec", np.linalg.eigh(np.cov(Z.T))[1][:, -1]), ("single neuron 0", np.eye(n)[0]), ("single neuron 7", np.eye(n)[7])]:
    y = Z @ proj; y = (y - y.mean())/y.std()
    phi = [abs(np.mean(np.exp(1j*t*y))) for t in ts]
    print(f"{name:22s} |phi(t)|: " + " ".join(f"{p:.4f}" for p in phi))
# Walsh concentration of the final-layer code law: largest |E prod_{i in S} sign z_i| over random small sets S
S_ = np.sign(Z)
vals = []
for k in [1, 2, 3, 4]:
    v = [abs(np.mean(np.prod(S_[:, rng.choice(n, k, replace=False)], 1))) for _ in range(2000)]
    vals.append((k, np.mean(v), np.max(v)))
print("Walsh coefficients of the final code law (|E prod sign|) over 2000 random sets of size k: " + "; ".join(f"k={k}: mean {a:.3f}, max {b:.3f}" for k, a, b in vals))
