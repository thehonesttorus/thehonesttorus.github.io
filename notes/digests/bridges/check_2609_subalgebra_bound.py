# Numerical check (own): does Cor. III.5 / Lemma III.4 of arXiv:2609.38007 hold when the
# tensor factor B(H_R) (x) 1_E is replaced by a general finite-dim subalgebra
# N = (+)_i M_{k_i} (x) 1_{m_i}, with E_N its trace-preserving conditional expectation?
import numpy as np
from numpy.linalg import eigh
def expm(H):
    H = (H + H.conj().T)/2; w, U = eigh(H); return (U*np.exp(w))@U.conj().T
def logm(A):
    A = (A + A.conj().T)/2; w, U = eigh(A); return (U*np.log(w))@U.conj().T
rng = np.random.default_rng(1)

def rand_state(n, spread=2.0):
    G = rng.normal(size=(n,n)) + 1j*rng.normal(size=(n,n))
    H = (G + G.conj().T)/2
    r = expm(-spread*H/np.linalg.norm(H,2)); return r/np.trace(r).real

def fpow(A, p):
    w, U = eigh(A); return (U*w**p)@U.conj().T
def fexpit(A, t):
    w, U = eigh(A); return (U*np.exp(1j*t*np.log(w)))@U.conj().T

def make_EN(blocks):
    # blocks: list of (k, m); N = (+) M_k (x) 1_m ; E_N = (+) id (x) tr_m/m on the diagonal blocks
    n = sum(k*m for k,m in blocks)
    def EN(X):
        Y = np.zeros_like(X); o = 0
        for k,m in blocks:
            B = X[o:o+k*m, o:o+k*m].reshape(k,m,k,m)
            red = np.einsum('ajbj->ab', B)/m
            Y[o:o+k*m, o:o+k*m] = np.kron(red, np.eye(m)); o += k*m
        return Y
    return n, EN

def relent(r, s): return np.trace(r@(logm(r)-logm(s))).real

def check(blocks, spread=2.0, perturb=None):
    n, EN = make_EN(blocks)
    rho = rand_state(n, spread)
    sig = rand_state(n, spread) if perturb is None else perturb(rho, EN)
    delta = relent(rho, sig) - relent(EN(rho), EN(sig))
    ts = np.linspace(-12, 12, 4801); ts = ts[np.abs(ts) > 1e-9]; dt = ts[1]-ts[0]
    eta = np.array([np.linalg.norm((lambda u: u-EN(u))(fexpit(sig,t)@fexpit(rho,-t)), 2) for t in ts])
    q = 0.5*np.sum(eta/np.abs(np.sinh(np.pi*ts)))*dt
    out = []
    for a in [0.25, 0.5, 1.0]:
        M = np.trace(fpow(rho, 1+a)@fpow(sig, -a)).real
        rhs = (1/a+3)*M**(1/(1+a))*q**(2*a/(1+a))
        out.append((a, rhs))
    return delta, q, out

for blocks in [[(2,1),(3,1)], [(2,2),(1,3)], [(3,2)], [(1,1),(1,1),(2,2)]]:
    for trial in range(3):
        d, q, out = check(blocks)
        ok = all(d <= rhs + 1e-9 for _, rhs in out)
        print(blocks, f"delta={d:.4e} q={q:.4e}", " ".join(f"a={a}:rhs={r:.3e}" for a,r in out), "OK" if ok else "VIOLATION")

# near-sufficient case: sigma's cocycle almost in N  (sigma = rho conjugated by small perturbation)
def near(rho, EN, eps=1e-2):
    n = rho.shape[0]
    G = rng.normal(size=(n,n)); G = (G+G.T)/2
    s = expm(logm(rho) + eps*G); return s/np.trace(s).real
for blocks in [[(2,2),(1,3)], [(3,2)]]:
    d, q, out = check(blocks, perturb=near)
    print("near", blocks, f"delta={d:.3e} q={q:.3e}", " ".join(f"a={a}:rhs={r:.3e}" for a,r in out))
