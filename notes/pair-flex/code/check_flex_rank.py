# The kernel of the symbols -> pair-slice map (Theorem 2) has dimension exactly C(K+3,4) + q(q-1)/2, q = K(K+1)/2,
# for generic loadings with n >= q (V held fixed). Builds the map's matrix on a basis and computes its nullity.
import itertools, math, numpy as np
rng = np.random.default_rng(3)
def sym4_basis(K):
    out = []
    for c in itertools.combinations_with_replacement(range(K), 4):
        Z = np.zeros((K,) * 4)
        for p in set(itertools.permutations(c)):
            Z[p] = 1.0
        out.append(Z)
    return out
def sym2_basis(K):
    out = []
    for i, j in itertools.combinations_with_replacement(range(K), 2):
        Z = np.zeros((K, K)); Z[i, j] = Z[j, i] = 1.0; out.append(Z)
    return out
for K, n in ((2, 6), (2, 10), (3, 9), (3, 14), (4, 12), (4, 24)):
    q = K * (K + 1) // 2
    U = rng.normal(size=(n, K))
    pairs = [(a, b) for a in range(n) for b in range(a + 1, n)]
    ordered = [(a, b) for a in range(n) for b in range(n) if a != b]
    def slices(T4, N, G, X):
        t4 = lambda a, b, c, d: np.einsum("pqrs,p,q,r,s->", T4, U[a], U[b], U[c], U[d])
        diag = [4 * G[a] @ U[a] - 6 * U[a] @ N[a] @ U[a] + 3 * t4(a, a, a, a) + X[a] for a in range(n)]
        k22 = [U[b] @ N[a] @ U[b] + U[a] @ N[b] @ U[a] - t4(a, a, b, b) for a, b in pairs]
        k31 = [G[a] @ U[b] for a, b in ordered]
        return np.array(diag + k22 + k31)
    cols = []
    z4, zN, zG, zX = np.zeros((K,) * 4), np.zeros((n, K, K)), np.zeros((n, K)), np.zeros(n)
    for B in sym4_basis(K):
        cols.append(slices(B, zN, zG, zX))
    for a in range(n):
        for B in sym2_basis(K):
            N = zN.copy(); N[a] = B; cols.append(slices(z4, N, zG, zX))
    for a in range(n):
        for k in range(K):
            G = zG.copy(); G[a, k] = 1.0; cols.append(slices(z4, zN, G, zX))
    for a in range(n):
        X = zX.copy(); X[a] = 1.0; cols.append(slices(z4, zN, zG, X))
    M = np.array(cols).T
    s = np.linalg.svd(M, compute_uv=False); rank = int(np.sum(s > 1e-9 * s[0])); gap = s[rank - 1] / s[0], (s[rank] / s[0] if rank < len(s) else 0.0)
    pred = math.comb(K + 3, 4) + q * (q - 1) // 2
    print(f"K={K} n={n}: parameters {M.shape[1]}, pair data {M.shape[0]}, nullity {M.shape[1] - rank}, "
          f"predicted C(K+3,4) + q(q-1)/2 = {math.comb(K + 3, 4)} + {q * (q - 1) // 2} = {pred}; n(q+1) - q(q-1)/2 - n(n+1)/2 = {n * (q + 1) - q * (q - 1) // 2 - n * (n + 1) // 2}; sv gap {gap[0]:.1e} | {gap[1]:.1e}")
