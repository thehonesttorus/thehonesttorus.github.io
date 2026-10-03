# Neural Gabber-Galil lemma, walk form: lump the interleaved GG walk on depth-2n source cylinders onto
# network-induced classes, then Metropolize toward the gain-weighted class law. Compare gaps with the bound.
import numpy as np, itertools
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist
rng = np.random.default_rng(5)
n = 5; M = 2**n; N = 2*n
maps = [lambda x,y: (x, (y+2*x)%M), lambda x,y: (x, (y-2*x)%M), lambda x,y: (x, (y+2*x+1)%M), lambda x,y: (x, (y-2*x-1)%M),
        lambda x,y: ((x+2*y)%M, y), lambda x,y: ((x-2*y)%M, y), lambda x,y: ((x+2*y+1)%M, y), lambda x,y: ((x-2*y-1)%M, y)]
V = M*M
P = np.zeros((V, V))
for x in range(M):
    for y in range(M):
        for f in maps:
            u, v = f(x, y); P[x*M+y, u*M+v] += 1/8
gap_src = 1 - np.sort(np.linalg.eigvalsh((P+P.T)/2))[-2]
# source cylinders -> Cantor-structured inputs (digits interleaved back to a binary word) -> random He ReLU network
width, L, r = 256, 8, 0.6
G = rng.standard_normal((N, 2, width))
def word(x, y):
    b = []
    for k in range(n): b += [(x >> k) & 1, (y >> k) & 1]
    return b
X = np.array([sum(r**k * G[k, b] for k, b in enumerate(word(x, y))) for x in range(M) for y in range(M)])
X /= np.linalg.norm(X, axis=1, keepdims=True)
H = X
for _ in range(L): H = np.maximum(H @ (rng.standard_normal((width, width)) * np.sqrt(2/width)).T, 0)
gain = np.abs(H).sum(1)                         # S = ||h_L||_1, the network gain on each source cylinder
Hn = H / np.linalg.norm(H, axis=1, keepdims=True)
for eps in [0.02, 0.05, 0.1]:
    cls = fcluster(linkage(pdist(Hn), method='single'), eps, criterion='distance') - 1   # PB eps-chain classes of the network output
    K = cls.max() + 1
    Pi = np.zeros((V, K)); Pi[np.arange(V), cls] = 1
    mu = np.full(V, 1/V); mub = Pi.T @ mu
    Pbar = (Pi.T * mu) @ P @ Pi / mub[:, None]                       # lumped chain (source stationary law)
    Pbar_s = np.diag(np.sqrt(mub)) @ Pbar @ np.diag(1/np.sqrt(mub))
    gap_lump = 1 - np.sort(np.linalg.eigvalsh((Pbar_s + Pbar_s.T)/2))[-2] if K > 1 else np.nan
    h = (Pi.T @ (mu*gain)) / mub; h = h / (mub @ h)                    # class-averaged normalized gain (target density)
    Kmet = Pbar * np.minimum(1, h[None, :] / h[:, None]); np.fill_diagonal(Kmet, 0); np.fill_diagonal(Kmet, 1 - Kmet.sum(1))
    nub = mub * h
    Ks = np.diag(np.sqrt(nub)) @ Kmet @ np.diag(1/np.sqrt(nub))
    gap_met = 1 - np.sort(np.linalg.eigvalsh((Ks + Ks.T)/2))[-2] if K > 1 else np.nan
    print(f"eps={eps}: {K:4d} network classes | gap source={gap_src:.4f}  lumped={gap_lump:.4f} (>= source)  "
          f"Metropolized to gain law={gap_met:.4f}  certified bound={gap_src*h.min()/h.max():.4f} (h spread {h.max()/h.min():.2f})")
