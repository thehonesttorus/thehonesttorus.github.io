# brute-force check of mclegs' all-distinct contractions on the width-7 net
import numpy as np, math
Wcol = np.load("../official/W_off9.npy").astype(np.float32); L, n, _ = Wcol.shape
F = np.load("mt_off9_full.npz"); G = np.load("mt_off9_legs.npz")
mu = F["mu"].astype(np.float64); sd = np.sqrt(F["var"].astype(np.float64))
Phi = 0.5 * (1 + np.vectorize(math.erf)(mu / sd / math.sqrt(2)))
N = 400000; rng = np.random.default_rng(11); bs = []; done = 0
while done < N:
    h = 0 if done < N // 2 else 1
    b = min(32768, N - done, (N // 2 - done) if h == 0 else N - done); bs.append(rng.standard_normal((n, b), dtype=np.float32)); done += b
Y = np.concatenate(bs, 1)
for l in range(L - 1):
    Z = Wcol[l] @ Y; Y = np.maximum(Z, 0)
    d = Z.astype(np.float64) - mu[l][:, None]
    T = np.einsum("ap,bp,cp->abc", d, d, d) / d.shape[1]
    idx = np.indices((n, n, n)); dist = (idx[0] != idx[1]) & (idx[1] != idx[2]) & (idx[0] != idx[2]); Td = np.where(dist, T, 0)
    X = Wcol[l + 1].astype(np.float64) * Phi[l][None, :]
    bf_d3 = np.einsum("ia,ib,ic,abc->i", X, X, X, Td)
    bf_d21 = np.einsum("ia,ib,cd,abd->ic", X, X, X, Td); np.fill_diagonal(bf_d21, 0)
    bf_allD = np.einsum("ib,ic,abc->ai", X, X, Td)
    XX = X * X; k3z = F["k3"][l].astype(np.float64); D21z = F["D21"][l].astype(np.float64).copy(); np.fill_diagonal(D21z, 0)
    u3 = G["u3"][l]; S = G["S"][l].astype(np.float64); P = G["P"][l].astype(np.float64)
    tri_d3 = u3 - 3 * np.einsum("ia,ai->i", XX, D21z @ X.T) - (XX * X) @ k3z
    tri_d21 = P - XX @ D21z @ X.T - 2 * (X * (D21z @ X.T).T) @ X.T - (XX * k3z[None, :]) @ X.T; np.fill_diagonal(tri_d21, 0)
    allD = S - (XX @ D21z).T - 2 * X.T * (D21z @ X.T) - (XX * k3z[None, :]).T
    r = lambda a, b: np.abs(a - b).max() / np.abs(b).max()
    print(l, "d3 %.1e d21 %.1e allD %.1e" % (r(tri_d3, bf_d3), r(tri_d21, bf_d21), r(allD, bf_allD)))
