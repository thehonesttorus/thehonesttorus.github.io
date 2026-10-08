# Exact check (finite enumeration, no sampling) of the conditional-independence identities and the flex theorem.
# Model: t takes finitely many values in R^K (arbitrary, non-Gaussian weights); given t, units a = 0..n-1 are independent,
# e_a = u_a . t + eps_a with eps_a | t a 3-point law whose variance/skewness depend on t (heteroscedastic, skewed).
import itertools, numpy as np
rng = np.random.default_rng(7)
n, K = 5, 2
tv = rng.normal(size=(4, K)); pt = rng.dirichlet(np.ones(4))          # law of t: 4 atoms
U = rng.normal(size=(n, K))
# per unit and per t-atom: eps support {x1,x2,x3} with probs, mean 0
eps_sup = np.zeros((n, 4, 3)); eps_p = np.zeros((n, 4, 3))
for a in range(n):
    for j in range(4):
        x = rng.normal(size=3) * (1 + 0.5 * np.tanh(tv[j] @ rng.normal(size=K)))
        p = rng.dirichlet(np.ones(3)); x = x - p @ x
        eps_sup[a, j], eps_p[a, j] = x, p
# joint law of (t, e): enumerate atoms
states, probs = [], []
for j in range(4):
    for idx in itertools.product(range(3), repeat=n):
        e = U @ tv[j] + np.array([eps_sup[a, j, idx[a]] for a in range(n)])
        pr = pt[j] * np.prod([eps_p[a, j, idx[a]] for a in range(n)])
        states.append(np.concatenate([e, tv[j]])); probs.append(pr)
X = np.array(states); P = np.array(probs); X = X - P @ X              # centre e and t
def cum4(i, j, k, l):
    m = lambda *s: P @ np.prod(X[:, list(s)], axis=1)
    return m(i, j, k, l) - m(i, j) * m(k, l) - m(i, k) * m(j, l) - m(i, l) * m(j, k)
T = lambda k: n + k
err = []
for a, b, c, d in itertools.permutations(range(n), 4):
    lhs = cum4(a, b, c, d)
    rhs = sum(U[a, p] * U[b, q] * U[c, r] * U[d, s] * cum4(T(p), T(q), T(r), T(s)) for p, q, r, s in itertools.product(range(K), repeat=4))
    err.append(abs(lhs - rhs) / (abs(lhs) + 1e-15))
print("(1,1,1,1) = T4[u_a,u_b,u_c,u_d]          max rel err %.1e" % max(err))
err = []
for a, c, d in itertools.permutations(range(n), 3):
    lhs = cum4(a, a, c, d)
    N = np.array([[cum4(a, a, T(p), T(q)) for q in range(K)] for p in range(K)])
    err.append(abs(lhs - U[c] @ N @ U[d]) / abs(lhs))
print("(2,1,1) = u_c^T kappa(e_a,e_a,t,t) u_d    max rel err %.1e" % max(err))
err = []
for a, b in itertools.permutations(range(n), 2):
    g = np.array([cum4(a, a, a, T(p)) for p in range(K)])
    err.append(abs(cum4(a, a, a, b) - g @ U[b]) / abs(cum4(a, a, a, b)))
print("(3,1)  = kappa(e_a,e_a,e_a,t) . u_b       max rel err %.1e" % max(err))
# (2,2) is NOT a pure evaluation: residual = K22 - u_b N_a u_b - u_a N_b u_a + T4[u_a^2 u_b^2] = Cov(v_a(t), v_b(t))
a, b = 0, 1
Na = np.array([[cum4(a, a, T(p), T(q)) for q in range(K)] for p in range(K)])
Nb = np.array([[cum4(b, b, T(p), T(q)) for q in range(K)] for p in range(K)])
T4ab = sum(U[a, p] * U[a, q] * U[b, r] * U[b, s] * cum4(T(p), T(q), T(r), T(s)) for p, q, r, s in itertools.product(range(K), repeat=4))
va = np.array([eps_p[a, j] @ eps_sup[a, j] ** 2 for j in range(4)]); vb = np.array([eps_p[b, j] @ eps_sup[b, j] ** 2 for j in range(4)])
covv = pt @ (va * vb) - (pt @ va) * (pt @ vb)
print("(2,2) residual vs Cov(v_a, v_b)          %.6e vs %.6e" % (cum4(a, a, b, b) - U[b] @ Na @ U[b] - U[a] @ Nb @ U[a] + T4ab, covv))
# diagonal: kappa4(e_a) = 4 g_a.u_a - 6 u_a^T N_a u_a + 3 T4[u_a^4] + 3 Var(v_a) + E kappa4(eps_a | t)
err = []
for a in range(n):
    Na = np.array([[cum4(a, a, T(p), T(q)) for q in range(K)] for p in range(K)])
    g = np.array([cum4(a, a, a, T(p)) for p in range(K)])
    T4aa = sum(U[a, p] * U[a, q] * U[a, r] * U[a, s] * cum4(T(p), T(q), T(r), T(s)) for p, q, r, s in itertools.product(range(K), repeat=4))
    v = np.array([eps_p[a, j] @ eps_sup[a, j] ** 2 for j in range(4)])
    x = np.array([eps_p[a, j] @ eps_sup[a, j] ** 4 - 3 * (eps_p[a, j] @ eps_sup[a, j] ** 2) ** 2 for j in range(4)])
    rhs = 4 * g @ U[a] - 6 * U[a] @ Na @ U[a] + 3 * T4aa + 3 * (pt @ (v * v) - (pt @ v) ** 2) + pt @ x
    err.append(abs(cum4(a, a, a, a) - rhs) / abs(cum4(a, a, a, a)))
print("diag   = 4 g.u - 6 uNu + 3 T4 + 3 Var v + E x   max rel err %.1e" % max(err))
# flex: a random dT4 in Sym^4(R^K) with dN_a = dT4[u_a,u_a,.,.]/2, and a random antisymmetric A on Sym^2(R^K);
# the pair-slice formulas are unchanged, the omitted classes are not.
from itertools import permutations
def symK4(Z):
    return sum(np.transpose(Z, p) for p in permutations(range(4))) / 24
dT = symK4(rng.normal(size=(K,) * 4))
iu = np.triu_indices(K); q = len(iu[0]); A = rng.normal(size=(q, q)); A = A - A.T
def vecS(S):                       # Frobenius-isometric coordinates on Sym^2
    w = np.where(iu[0] == iu[1], 1.0, np.sqrt(2.0)); return S[iu] * w
def matS(v):
    w = np.where(iu[0] == iu[1], 1.0, 1 / np.sqrt(2.0)); S = np.zeros((K, K)); S[iu] = v * w; return S + np.triu(S, 1).T
dN1 = np.array([0.5 * np.einsum("pqrs,p,q->rs", dT, U[a], U[a]) for a in range(n)])
dN2 = np.array([matS(A @ vecS(np.outer(U[a], U[a]))) for a in range(n)])
T4f = lambda Z, x, y, z, w: np.einsum("pqrs,p,q,r,s->", Z, x, y, z, w)
for name, dN, dT4 in (("Sym^4 flex", dN1, dT), ("antisymmetric flex", dN2, 0 * dT)):
    dK22 = max(abs(U[b] @ dN[a] @ U[b] + U[a] @ dN[b] @ U[a] - T4f(dT4, U[a], U[a], U[b], U[b])) for a, b in permutations(range(n), 2))
    ddiag = max(abs(-6 * U[a] @ dN[a] @ U[a] + 3 * T4f(dT4, U[a], U[a], U[a], U[a])) for a in range(n))
    d211 = max(abs(U[c] @ dN[a] @ U[d]) for a, c, d in permutations(range(n), 3))
    d1111 = max(abs(T4f(dT4, U[a], U[b], U[c], U[d])) for a, b, c, d in permutations(range(n), 4))
    print(f"{name:20s}: change of K22 {dK22:.1e}, of diag {ddiag:.1e} (K31 untouched); change of (2,1,1) {d211:.2e}, of (1,1,1,1) {d1111:.2e}")
