import numpy as np, itertools
rng = np.random.default_rng(7); n = 6
W = rng.standard_normal((n, n)) / np.sqrt(n); Ph = rng.random(n); rho = rng.random(n); c13 = rng.standard_normal(n)
A = rng.standard_normal((n, n)); C = A + A.T; np.fill_diagonal(C, 0); Gm = rng.standard_normal((n, n)); np.fill_diagonal(Gm, 0)
def kap(a, b, c):
    s = 0.0
    for (u, v, w) in itertools.permutations((a, b, c)):
        s += Gm[u, v] / 2 * c13[u] * Ph[v] * Ph[w] * C[u, w]          # GC1
    return s
def kap2(a, b, c):
    return sum(Gm[u, v] / 2 * rho[u] * rho[v] * Ph[w] * C[v, w] for (u, v, w) in itertools.permutations((a, b, c)))
def kds(a, b, c):
    return sum(C[u, v] ** 2 / 2 * c13[u] * rho[v] * Ph[w] * C[u, w] for (u, v, w) in itertools.permutations((a, b, c)))
def ktri(a, b, c): return rho[a] * rho[b] * rho[c] * C[a, b] * C[b, c] * C[c, a]
bf = {k: np.zeros(n) for k in ("gc1", "gc2", "ds", "tri")}
for i in range(n):
    for a, b, c in itertools.permutations(range(n), 3):
        w = W[i, a] * W[i, b] * W[i, c]
        bf["gc1"][i] += w * kap(a, b, c); bf["gc2"][i] += w * kap2(a, b, c); bf["ds"][i] += w * kds(a, b, c); bf["tri"][i] += w * ktri(a, b, c)
X = W * Ph[None, :]; XX = X * X; Wr = W * rho[None, :]; CX = C @ X.T; GX = Gm @ X.T
gc1 = 3 * np.einsum("iu,ui->i", W * c13[None, :], GX * CX - (Gm * C) @ XX.T)
Y = Wr * CX.T
gc2 = 3 * (np.einsum("iu,ui->i", Wr, Gm @ Y.T) - np.einsum("iu,ui->i", Wr * X, (Gm * C) @ Wr.T))
Q = Wr @ (C * C)
ds = 3 * np.einsum("iu,ui->i", W * c13[None, :], Q.T * CX - (C * C * C) @ (Wr * X).T)
tri = np.array([Wr[i] @ ((C * ((C * Wr[i][None, :]) @ C)) @ Wr[i]) for i in range(n)])
for k, v in (("gc1", gc1), ("gc2", gc2), ("ds", ds), ("tri", tri)):
    print(k, "%.1e" % (np.abs(v - bf[k]).max() / np.abs(bf[k]).max()))
