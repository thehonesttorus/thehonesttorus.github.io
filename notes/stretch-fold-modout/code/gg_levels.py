# Gabber-Galil graphs on (Z/p^n)^2 = level-n cylinder quotients of the p-adic plane Z_p^2.
# (1) second eigenvalue uniform in n (local-to-global expansion at every level);
# (2) the tail operator M_n (GG maps applied to the coordinate below level n) fixes level-n functions and
#     has gap >= 1 - lambda* on their orthocomplement, so (1-lambda*)(I-E_n) <= I-M_n <= 2(I-E_n).
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
def gg_maps(m):
    return [lambda x,y: (x, (y+2*x)%m), lambda x,y: (x, (y-2*x)%m),
            lambda x,y: (x, (y+2*x+1)%m), lambda x,y: (x, (y-2*x-1)%m),
            lambda x,y: ((x+2*y)%m, y), lambda x,y: ((x-2*y)%m, y),
            lambda x,y: ((x+2*y+1)%m, y), lambda x,y: ((x-2*y-1)%m, y)]
def markov(m):
    X, Y = np.meshgrid(np.arange(m), np.arange(m), indexing='ij'); X = X.ravel(); Y = Y.ravel()
    src = X*m + Y; rows, cols = [], []
    for f in gg_maps(m):
        u, v = f(X, Y); rows.append(src); cols.append(u*m + v)
    P = sp.csr_matrix((np.full(8*m*m, 1/8), (np.concatenate(rows), np.concatenate(cols))), shape=(m*m, m*m))
    return P
lam_star = 5*np.sqrt(2)/8
print(f"Gabber-Galil bound lambda* = 5*sqrt(2)/8 = {lam_star:.4f}")
for p, ns in [(2, range(1, 8)), (3, range(1, 5)), (5, range(1, 4)), (7, range(1, 3))]:
    out = []
    for n in ns:
        m = p**n; P = markov(m); S = (P + P.T) / 2
        if m*m <= 400:
            ev = np.sort(np.linalg.eigvalsh(S.toarray()))
        else:
            ev = np.sort(sla.eigsh(S, k=3, which='LA', return_eigenvectors=False))
        out.append(f"{m}:{ev[-2]:.4f}")
    print(f"p={p}: second eigenvalue of the symmetrized walk on (Z/p^n)^2  ->  " + "  ".join(out))
# tail operator check on Z/p^N with level n: functions on (Z/p^N)^2, M_n acts on w where v = r + p^n w
p, N, n = 2, 6, 2
m = p**N; k = p**(N-n)
X, Y = np.meshgrid(np.arange(m), np.arange(m), indexing='ij'); X = X.ravel(); Y = Y.ravel()
rx, ry, wx, wy = X % p**n, Y % p**n, X // p**n, Y // p**n
rows, cols = [], []
for f in gg_maps(k):
    u, v = f(wx, wy); rows.append(X*m + Y); cols.append((rx + p**n*u)*m + (ry + p**n*v))
Mn = sp.csr_matrix((np.full(8*m*m, 1/8), (np.concatenate(rows), np.concatenate(cols))), shape=(m*m, m*m)).toarray()
Mn = (Mn + Mn.T) / 2
lev = (X % p**n) * p**n + (Y % p**n)                     # level-n cylinder label
En = (lev[:, None] == lev[None, :]).astype(float); En /= En.sum(1, keepdims=True)
print(f"tail operator check (p={p}, N={N}, level n={n}): ||M_n E_n - E_n|| = {np.abs(Mn@En - En).max():.1e}, ||[M_n,E_n]|| = {np.abs(Mn@En - En@Mn).max():.1e}")
ev = np.linalg.eigvalsh((np.eye(m*m) - En) @ Mn @ (np.eye(m*m) - En))
ev = ev[np.abs(ev) > 1e-9]
print(f"  spectrum of M_n on the orthocomplement of level-n functions: [{ev.min():.4f}, {ev.max():.4f}]  (gap >= 1 - lambda* = {1-lam_star:.4f})")
