# closed-form transport of the one-facet second-order kappa3 class against brute force (algebra only)
import numpy as np
rng = np.random.default_rng(4); n = 7
W = rng.standard_normal((n, n)) / np.sqrt(n); Ph = rng.random(n); rho = rng.random(n)
A = rng.standard_normal((n, n)); C = A + A.T; np.fill_diagonal(C, 0)
X = W * Ph[None, :]
def kap(a, b, c):
    return rho[a] * Ph[b] * Ph[c] * C[a, b] * C[a, c] + Ph[a] * rho[b] * Ph[c] * C[a, b] * C[b, c] + Ph[a] * Ph[b] * rho[c] * C[a, c] * C[b, c]
bf3 = np.zeros(n); bf21 = np.zeros((n, n))
for i in range(n):
    for a in range(n):
        for b in range(n):
            for d in range(n):
                if len({a, b, d}) < 3: continue
                k = kap(a, b, d); bf3[i] += W[i, a] * W[i, b] * W[i, d] * k
                for c in range(n): bf21[i, c] += W[i, a] * W[i, b] * W[c, d] * k
np.fill_diagonal(bf21, 0)
G = C @ X.T                                   # G_ai = sum_b C_ab X_ib
XX = X * X; C2 = C * C; Wr = W * rho[None, :]
f3 = 3 * np.einsum("ia,ia->i", Wr, G.T ** 2 - XX @ C2)
M1 = G * G - C2 @ XX.T                         # (d, i)
Q = Wr @ C2                                    # Q_ib = sum_a W_ia rho_a C_ab^2
f21 = (Wr @ M1).T + 2 * ((Wr * G.T) @ G - (X * Q) @ X.T)
np.fill_diagonal(f21, 0)
print("facet D3 %.1e  D21 %.1e" % (np.abs(f3 - bf3).max() / np.abs(bf3).max(), np.abs(f21 - bf21).max() / np.abs(bf21).max()))
