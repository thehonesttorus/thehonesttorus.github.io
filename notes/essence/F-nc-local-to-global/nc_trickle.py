"""Check of Theorem F1 (k-level spectral trickle-down in a forget lattice) on a non-commuting quantum Gibbs state.
Hilbert space: L2(M_2^{(x)n}) with the KMS inner product <X,Y> = Tr(X^* s^{1/2} Y s^{1/2}), s = exp(-bH)/Z,
H a random Heisenberg-plus-field ring (non-commuting). Forget lattice: N_A = operators acting trivially on A;
E_A = KMS-orthogonal projection onto N_A. Prints, per level m, eta_m = max_A ||sum_{i in A}(E_{A-i}-E_A)|| - 1,
the pairwise-defect bound lambda_max(C^A), the theorem's gap bound prod_m (1-(1+eta_m)/m), and the true gap of
(1/n) sum_i E_i on the orthocomplement of E_V."""
import itertools, numpy as np, sys
n = int(sys.argv[1]) if len(sys.argv) > 1 else 4; beta = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
rng = np.random.default_rng(1); d = 2 ** n
X = np.array([[0, 1], [1, 0]], complex); Y = np.array([[0, -1j], [1j, 0]]); Z = np.diag([1., -1]).astype(complex); I2 = np.eye(2)
def op(P, i):
    out = np.array([[1.]]);
    for k in range(n): out = np.kron(out, P if k == i else I2)
    return out
H = sum(rng.normal(1, .3) * op(P, i) @ op(P, (i + 1) % n) for i in range(n) for P in (X, Y, Z))
H = H + sum(rng.normal(0, .7) * op(P, i) for i in range(n) for P in (X, Z))
w, V = np.linalg.eigh(H); s = V @ np.diag(np.exp(-beta * (w - w.min()))) @ V.conj().T; s /= np.trace(s).real
w2, V2 = np.linalg.eigh(s); sq = V2 @ np.diag(np.sqrt(w2)) @ V2.conj().T
# KMS isometry: X -> s^{1/4} X s^{1/4} maps KMS inner product to Hilbert-Schmidt; use vec of sqrt-sandwich
w4 = V2 @ np.diag(w2 ** .25) @ V2.conj().T
def sandwich(Xop): return (w4 @ Xop @ w4).reshape(-1)
paulis = [I2, X, Y, Z]
def subspace(A):  # orthonormal basis (in HS after sandwich) of operators trivial on A
    cols = []
    for idx in itertools.product(range(4), repeat=n):
        if any(idx[a] != 0 for a in A): continue
        M = np.array([[1.]])
        for k in range(n): M = np.kron(M, paulis[idx[k]])
        cols.append(sandwich(M))
    Q, _ = np.linalg.qr(np.array(cols).T); return Q
proj = {}
def E(A):
    A = tuple(sorted(A))
    if A not in proj: Q = subspace(A); proj[A] = Q @ Q.conj().T
    return proj[A]
eta = {}
for m in range(2, n + 1):
    e_max = c_max = 0
    for A in itertools.combinations(range(n), m):
        Ps = [E(set(A) - {i}) - E(A) for i in A]
        e_max = max(e_max, np.linalg.norm(sum(Ps), 2) - 1)
        C = np.array([[0 if i == j else np.linalg.norm(Ps[i] @ Ps[j], 2) for j in range(m)] for i in range(m)])
        c_max = max(c_max, np.linalg.eigvalsh(C).max())
    eta[m] = e_max; print(f"m={m}: eta_m={e_max:.4f}  lambda_max(defects)={c_max:.4f}  factor 1-(1+eta)/m={1-(1+e_max)/m:.4f}")
bound = np.prod([1 - (1 + eta[m]) / m for m in range(2, n + 1)])
Pg = sum(E({i}) for i in range(n)) / n; EV = E(range(n)); Id = np.eye(d * d)
# gap of (I - Pg) restricted to range(I - EV)
Qc = np.linalg.svd(Id - EV)[0][:, : d * d - 1]
ev = np.linalg.eigvalsh(Qc.conj().T @ (Id - Pg) @ Qc)
print(f"n={n} beta={beta}: theorem gap bound {bound:.4f}   true gap {ev.min():.4f}   ratio {ev.min()/bound:.2f}")
