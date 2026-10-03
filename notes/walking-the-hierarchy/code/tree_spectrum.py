# Exact spectrum of last-point-of-contact (ultrametric radial) kernels on a random weighted tree:
#   K f(x) = sum_y psi(x ^ y) f(y) mu(y)  (x^y = last common vertex), Markov after normalisation.
# Claim: each detail space D_v (functions of v's children with mu-mean zero, supported in v) is an eigenspace with
#   lambda_v = 1 - Abar_v - psi(v) mu([v]),   Abar_v = sum_{u strict ancestor of v} psi(u) mu(shell_u toward v).
# Also: exact variance of walk averages (1/T) sum_v e_v g_T(lambda_v) vs simulation.
import numpy as np
rng = np.random.default_rng(1)
# random tree: depth 4, branching 2..4, random child masses
def build(depth):
    leaves = [((), 1.0)]
    for _ in range(depth):
        new = []
        for path, m in leaves:
            b = rng.integers(2, 5); p = rng.dirichlet(np.ones(b))
            new += [(path + (i,), m*p[i]) for i in range(b)]
        leaves = new
    return leaves
leaves = build(4); P = [l[0] for l in leaves]; mu = np.array([l[1] for l in leaves]); M = len(P)
def meet(a, b):
    k = 0
    while k < len(a) and a[k] == b[k]: k += 1
    return a[:k]
verts = sorted({p[:k] for p in P for k in range(len(p) + 1)}, key=len)
psi0 = {v: rng.uniform(0.2, 3.0) for v in verts if len(v) < 4}          # psi on internal vertices
# mass of a vertex
mass = {v: mu[[p[:len(v)] == v for p in P]].sum() for v in verts}
Kraw = np.array([[psi0.get(meet(P[i], P[j]), 0.0)*mu[j] if i != j else 0.0 for j in range(M)] for i in range(M)])
# Markov normalisation: add holding so that rows sum to 1 (scale psi by c so max row sum <= 1)
c = 1.0/Kraw.sum(1).max(); psi = {v: c*w for v, w in psi0.items()}
K = c*Kraw; K[np.diag_indices(M)] = 1 - K.sum(1)                         # holding = "contact at the leaf itself"
# predicted eigenvalue on D_v: 1 - Abar_v - psi(v) mass(v)  (holding counts as contact below v)
def lam(v):
    Ab = sum(psi[v[:k]]*(mass[v[:k]] - mass[v[:k+1]]) for k in range(len(v)))
    return 1 - Ab - psi[v]*mass[v]
# test: build D_v basis vectors and check K g = lambda_v g
errs = []
for v in verts:
    if len(v) == 4: continue
    kids = sorted({p[:len(v)+1] for p in P if p[:len(v)] == v})
    for a in range(len(kids) - 1):
        g = np.zeros(M)
        ia = [i for i, p in enumerate(P) if p[:len(v)+1] == kids[a]]; ib = [i for i, p in enumerate(P) if p[:len(v)+1] == kids[a+1]]
        g[ia] = 1/mass[kids[a]]; g[ib] = -1/mass[kids[a+1]]               # mu-mean zero over v's children
        errs.append(np.abs(K @ g - lam(v)*g).max())
print(f"tree with {M} leaves, {len(verts)} vertices: max |K g - lambda_v g| over all detail vectors = {max(errs):.2e}")
# reversibility
print(f"detailed balance max |mu_i K_ij - mu_j K_ji| = {np.abs(mu[:,None]*K - (mu[:,None]*K).T).max():.2e}")
# walk variance: exact formula vs simulation for f random
f = rng.standard_normal(M); T = 50
ev, V = np.linalg.eig(K.T)                                                # check via spectral sum on mu-weighted space
Dh = np.sqrt(mu); S = Dh[:, None]*K/Dh[None, :]                           # symmetric similar matrix
w, U = np.linalg.eigh((S + S.T)/2)
fc = (f - (mu*f).sum())*Dh; coef = U.T @ fc
def gT(l):
    if abs(1 - l) < 1e-12: return T
    return (1 + l)/(1 - l) - 2*l*(1 - l**T)/(T*(1 - l)**2)
pred = sum(coef[k]**2*gT(w[k]) for k in range(M))/T
R = 40000; cum = np.cumsum(K, 1)
x = rng.choice(M, size=R, p=mu); acc = np.zeros(R)
for t in range(T):
    acc += f[x]; u = rng.random(R); x = (u[:, None] > cum[x]).sum(1)
emp = (acc/T).var()
print(f"walk-average variance (T={T}): exact spectral formula {pred:.5f} vs simulation {emp:.5f};  iid would give {((f-(mu*f).sum())**2*mu).sum()/T:.5f}")
