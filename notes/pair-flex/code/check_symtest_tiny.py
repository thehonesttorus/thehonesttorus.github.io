# Bookkeeping check for symtest.py on a width-7 network, same samples, float64:
#  (a) R_true = W#[omitted classes of the empirical kappa4(y)] equals t' - T_pair(y slices);
#  (b) symbols() reproduces the brute-force projections U^T Psi_a^off U and U^(x4) Omega at K = 4 and K = 7;
#  (c) at K = n, R_pred = R_true exactly.
import sys, itertools, numpy as np
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from symtest import symbols, T_pair, r_pred
rng = np.random.default_rng(11)
n, N, l = 7, 400000, 1
Wcol = rng.normal(size=(3, n, n)) * np.sqrt(2.0 / n)
X = rng.normal(size=(n, N)); Z0 = Wcol[0] @ X; Y0 = np.maximum(Z0, 0); Z1 = Wcol[1] @ Y0; Y1 = np.maximum(Z1, 0)
W = Wcol[2]; H = W * W; Zp = W @ Y1
y = Y1; ybar = y.mean(1)
def kappa4_tensor(V):
    x = V - V.mean(1, keepdims=True); M = len(x)
    C = x @ x.T / x.shape[1]
    m4 = np.einsum("an,bn,cn,dn->abcd", x, x, x, x, optimize=True) / x.shape[1]
    return m4 - (np.einsum("ab,cd->abcd", C, C) + np.einsum("ac,bd->abcd", C, C) + np.einsum("ad,bc->abcd", C, C)), C
k4y, Cy = kappa4_tensor(y)
nd = np.zeros((n,) * 4, int)
for idx in itertools.product(range(n), repeat=4):
    nd[idx] = len(set(idx))
omit = np.where(nd >= 3, k4y, 0.0)
Tz = np.einsum("ia,jb,kc,ld,abcd->ijkl", W, W, W, W, omit, optimize=True)
ia, ib = np.array([a for a in range(n) for b in range(n) if a != b]), np.array([b for a in range(n) for b in range(n) if a != b])
def sl(T):
    k31 = np.array([[T[i, i, i, j] for j in range(n)] for i in range(n)]); np.fill_diagonal(k31, 0.0)
    return np.array([T[i, i, i, i] for i in range(n)]), np.array([T[i, i, j, j] for i, j in zip(ia, ib)]), k31
Rtrue = sl(Tz)
k4z, _ = kappa4_tensor(Zp)
t = sl(k4z)
K22y = np.array([[k4y[a, a, b, b] for b in range(n)] for a in range(n)]); np.fill_diagonal(K22y, 0.0)
K31y = np.array([[k4y[a, a, a, b] for b in range(n)] for a in range(n)]); np.fill_diagonal(K31y, 0.0)
d4y = np.array([k4y[a, a, a, a] for a in range(n)])
R2 = tuple(a - b for a, b in zip(t, T_pair(W, H, ia, ib, d4y, K22y, K31y)))
print("(a) t' - T_pair(y) vs W#[omitted]:", " ".join(f"{np.linalg.norm(a - b) / np.linalg.norm(a):.1e}" for a, b in zip(Rtrue, R2)))
# mcsym-style sums around a shifted centre (exercises the mean corrections)
ev, V = np.linalg.eigh(Cy); Ufull = V[:, ::-1]
Kmax = n; iu = np.triu_indices(Kmax)
mhat = ybar + 0.013
e = y - mhat[:, None]; E2 = e * e; E3 = E2 * e; tt = Ufull.T @ e; P = tt[iu[0]] * tt[iu[1]]
UU = (Ufull[:, iu[0]] * Ufull[:, iu[1]]).T; S = UU @ E2
acc = {"c": np.array(N), "K": np.array(Kmax), "layers": np.array([l]), f"U_{l}": Ufull}
for k, v in dict(s1=e.sum(1), s2=E2.sum(1), s3=E3.sum(1), s4=(E2 * E2).sum(1), C=e @ e.T, D=E2 @ e.T, et=e @ tt.T, e2t=E2 @ tt.T,
                 e3t=E3 @ tt.T, etp=e @ P.T, e2tp=E2 @ P.T, e2S=E2 @ S.T, t1=tt.sum(1), tp=P.sum(1), tpt=P @ tt.T, tptp=P @ P.T,
                 S1=S.sum(1), SS=S @ S.T).items():
    acc[f"{k}_{l}"] = v
for K in (4, 7):
    Sd = symbols(acc, l, K); U = Sd["U"]
    Noff_b = np.zeros((n, K, K))
    for a in range(n):
        Psi = np.array([[k4y[a, a, c, d] if (c != d and c != a and d != a) else 0.0 for d in range(n)] for c in range(n)])
        Noff_b[a] = U.T @ Psi @ U
    Om = np.where(nd == 4, k4y, 0.0)
    T4_b = np.einsum("abcd,ak,bm,cp,dq->kmpq", Om, U, U, U, U, optimize=True)
    print(f"(b) K={K}: N^off rel err {np.linalg.norm(Sd['Noff'] - Noff_b) / np.linalg.norm(Noff_b):.1e}, "
          f"T4^off rel err {np.linalg.norm(Sd['T4off'] - T4_b) / np.linalg.norm(T4_b):.1e}")
    Pd = r_pred(W, H, ia, ib, Sd)[0]
    print(f"    R_pred vs R_true (exact only at K = n): " + " ".join(f"{np.linalg.norm(a - b) / np.linalg.norm(a):.1e}" for a, b in zip(Rtrue, Pd)))
