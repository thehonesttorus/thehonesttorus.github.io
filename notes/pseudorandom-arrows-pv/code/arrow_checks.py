# Checks for the arrow-algebra theorems on A_N = M_2^{(x)N}, N small.
import numpy as np, itertools
rng = np.random.default_rng(7)
N = 5; d = 2**N
X = np.array([[0,1],[1,0]]); Z = np.diag([1,-1]); I2 = np.eye(2)
def weyl(a, b):
    M = np.array([[1.0]])
    for ai, bi in zip(a, b):
        M = np.kron(M, np.linalg.matrix_power(X, ai) @ np.linalg.matrix_power(Z, bi))
    return M
labels = list(itertools.product([0,1], repeat=2*N))
def depth(w):
    a, b = w[:N], w[N:]
    idx = [i+1 for i in range(N) if a[i] or b[i]]
    return max(idx) if idx else 0
lam = [0.0] + [(17 * 9**(m-1) - 1) / 4 for m in range(1, N+1)]      # canonical PB spectrum 4, 38, 344, ...
tau = lambda A: np.trace(A).real / d
def E(n, A):                                    # normalized partial trace onto the first n qubits
    T = A.reshape([2]*(2*N))
    k = N - n
    if k == 0: return A
    T = A.reshape(2**n, 2**k, 2**n, 2**k)
    R = np.einsum('ajbj->ab', T) / 2**k
    return np.kron(R, np.eye(2**k))
def T_heat(t, A):
    out = (1 - np.exp(-t*lam[1])) * E(0, A)
    for n in range(1, N+1):
        pn = np.exp(-t*lam[n]) - (np.exp(-t*lam[n+1]) if n < N else 0)
        out = out + pn * E(n, A)
    return out
# (1) Weyl diagonalization of the heat channel
t = 0.003
errs = []
for w in labels[::7]:
    W = weyl(w[:N], w[N:])
    errs.append(np.abs(T_heat(t, W) - np.exp(-t*lam[depth(w)]) * W).max())
print(f"(1) max |T_t W(w) - e^(-t lam_depth) W(w)| over sampled labels = {max(errs):.1e}")
# (2) sparse hierarchical Levy generator from random small-bias sets on the tail phase space
def symp(g, w):
    return (np.dot(g[:N], w[N:]) + np.dot(g[N:], w[:N])) % 2
def tail_set(n, k):
    S = []
    for _ in range(k):
        a = np.zeros(N, int); b = np.zeros(N, int)
        a[n:] = rng.integers(0, 2, N-n); b[n:] = rng.integers(0, 2, N-n)
        S.append(np.concatenate([a, b]))
    return S
k = 40
Sets = [tail_set(n, k) for n in range(N)]
eps = []
for n in range(N):
    biases = []
    for w in labels:
        w = np.array(w)
        if depth(w) > n:
            tailpart = np.concatenate([w[n:N], w[N+n:]])
            if tailpart.any():
                biases.append(abs(np.mean([(-1)**symp(g, w) for g in Sets[n]])))
    eps.append(max(biases))
print(f"(2) empirical bias of the random tail sets (k={k}) per level: {np.round(eps, 3)}")
psi_ratio = []
for w in labels:
    w = np.array(w); m = depth(w)
    if m == 0: continue
    psi = sum((lam[n+1]-lam[n]) * (1 - np.mean([(-1)**symp(g, w) for g in Sets[n]])) for n in range(N))
    psi_ratio.append(psi / lam[m])
e = max(eps)
print(f"    psi_hat/lam_depth in [{min(psi_ratio):.3f}, {max(psi_ratio):.3f}]  vs theorem bounds [1-eps, 1+eps] = [{1-e:.3f}, {1+e:.3f}]")
# build the sparse Lindbladian as a superoperator and check CP-semigroup facts on a random matrix
def Lhat(A):
    out = np.zeros_like(A, dtype=complex)
    for n in range(N):
        for g in Sets[n]:
            W = weyl(g[:N], g[N:]); out += (lam[n+1]-lam[n]) / k * (W @ A @ W.conj().T - A)
    return out
basis = [np.eye(d)[:, [i]] @ np.eye(d)[[j], :] for i in range(d) for j in range(d)]
Lmat = np.array([Lhat(B).reshape(-1) for B in basis]).T
from scipy.linalg import expm
A = rng.standard_normal((d, d)) + 1j * rng.standard_normal((d, d))
for t in [0.0005, 0.003, 0.02]:
    Th = (expm(t * Lmat) @ A.reshape(-1)).reshape(d, d)
    Tt = T_heat(t, A)
    rel = np.linalg.norm(Th - Tt) / np.linalg.norm(A - tau(A) * np.eye(d))
    print(f"    t={t}: ||T_hat(A)-T(A)||_2/||A-tau(A)||_2 = {rel:.3e}  (uniform bound eps/(e(1-eps)) = {e/(np.e*(1-e)):.3e});  "
          f"trace kept: {abs(np.trace(Th)-np.trace(A)):.1e}")
# Choi matrix positivity of exp(tL) for one t
t = 0.003; Phi = expm(t * Lmat)
choi = sum(np.kron(basis[i], (Phi @ basis[i].reshape(-1)).reshape(d, d)) for i in range(d*d))
print(f"    Choi(exp(tL_hat)) min eigenvalue = {np.linalg.eigvalsh((choi+choi.conj().T)/2).min():.2e} (>= 0: CP)")
# (3) renormalization identity  T_t o iota = iota o [e^{-2t} T_{9t} + (1-e^{-2t}) E_0]  (iota = 1 (x) . shift)
N0 = N - 1
B = rng.standard_normal((2**N0, 2**N0))
iotaB = np.kron(I2, B)
lhs = T_heat(t, iotaB)
def T_heat_small(t, A):                      # heat on N0 qubits with the same spectrum
    out = (1 - np.exp(-t*lam[1])) * np.trace(A)/A.shape[0] * np.eye(A.shape[0])
    for n in range(1, N0+1):
        pn = np.exp(-t*lam[n]) - (np.exp(-t*lam[n+1]) if n < N0 else 0)
        k2 = N0 - n
        R = np.einsum('ajbj->ab', A.reshape(2**n, 2**k2, 2**n, 2**k2)) / 2**k2
        out = out + pn * np.kron(R, np.eye(2**k2))
    return out
# note: on N qubits the shifted operator sits in qubits 2..N; truncation N vs N0 consistent since lam_{m+1}=9 lam_m+2
rhs = np.kron(I2, np.exp(-2*t) * T_heat_small(9*t, B) + (1-np.exp(-2*t)) * np.trace(B)/B.shape[0] * np.eye(B.shape[0]))
print(f"(3) renormalization identity residual = {np.abs(lhs - rhs).max():.1e}")
