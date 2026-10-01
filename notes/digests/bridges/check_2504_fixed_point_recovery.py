# Numerical check (own) for the digest arxiv-2504.02208.md, section 6.
#
# Claims checked, for a KMS-detailed-balanced generator L_A whose jumps are the single-site
# Paulis on A (here a Davies generator with Metropolis rates, which is KMS- and GNS-symmetric;
# Chen-Rouze use the quasi-local CKG23 generator instead, but the claims below only use
# KMS symmetry and "the Bohr components of the jumps sum to the Paulis on A"):
#   (a) Pauli twirl identity  rho - rho_{-A} = 2^{-2|A|-1} sum_{S in P_A} [S,[S,rho]]           (CR, proof of Thm III.1)
#   (b) Pimsner-Popa / index bound  rho <= 4^{|A|} rho_{-A}                                     (own, 1 line)
#   (c) KMS symmetry of L_A^dagger                                                              (CR Thm II.1 analogue)
#   (d) F_A := ker L_A^dagger  is contained in  N_A := 1_A (x) B(H_{A^c})                      (CR eq. (4.1))
#   (e) the t -> infinity limit P_A of the time average R_{A,t} recovers EXACTLY: P_A[rho_{-A}] = rho   (own derivation)
#   (f) F_A is contained in N_A^sigma := largest ad_H-invariant (= modular-invariant) subspace of N_A    (Takesaki, own use)
#   (g) commuting H: F_A large and P_A acts as the identity on operators far from A (local);
#       non-commuting H: F_A = C 1 here, so P_A is the global replacement channel (non-local)
#   (h) ||R_{A,t}[rho_{-A}] - rho||_1 -> 0 as t grows (in finite dimension at rate ~ 1/(gap t))
#   (i) CMI as a Petz sufficiency defect: I(A:C|B)_rho = D(rho||rho_{-A}) - D(rho_AB||(rho_{-A})_AB)   (own, 2 lines)
#   (j) non-commuting case: the extra fixed point is the ring reflection fixing the site in A (a symmetry acting off A)
import numpy as np
np.set_printoptions(precision=3, suppress=True)
X = np.array([[0, 1], [1, 0]], complex); Y = np.array([[0, -1j], [1j, 0]]); Z = np.diag([1., -1]).astype(complex); I2 = np.eye(2)

def op(single, site, n):
    out = np.array([[1.]], complex)
    for k in range(n):
        out = np.kron(out, single if k == site else I2)
    return out

def ham(n, J=1.0, g=0.0, h=0.3):
    H = np.zeros((2**n, 2**n), complex)
    for i in range(n):
        H += J*op(Z, i, n)@op(Z, (i+1) % n, n) + g*op(X, i, n) + h*op(Z, i, n)
    return H

def davies_heis(H, jumps, beta, tol=1e-9):
    """Heisenberg generator L^dag as a matrix on row-major vec, built in the computational basis."""
    E, U = np.linalg.eigh(H); d = len(E)
    Ldag = np.zeros((d*d, d*d), complex)
    nu = E[:, None] - E[None, :]                      # nu[i,j] = E_i - E_j : entry (i,j) of A_nu moves E_j -> E_i
    rate = np.minimum(1.0, np.exp(-beta*nu))          # Metropolis: raising energy by nu>0 costs e^{-beta nu}
    for A in jumps:
        At = U.conj().T@A@U                           # eigenbasis
        # group entries by Bohr frequency
        keys = np.round(nu/tol).astype(np.int64)
        for k in np.unique(keys):
            mask = (keys == k)
            Anu = np.where(mask, At, 0)
            if np.abs(Anu).max() < 1e-14: continue
            g = rate[mask][0]
            Anu_c = U@Anu@U.conj().T                  # back to computational basis
            AdA = Anu_c.conj().T@Anu_c
            # X -> g (Anu^dag X Anu - 1/2 {Anu^dag Anu, X});  vec(A X B) = kron(A, B.T) vec(X)
            Ldag += g*(np.kron(Anu_c.conj().T, Anu_c.T) - 0.5*np.kron(AdA, np.eye(d)) - 0.5*np.kron(np.eye(d), AdA.T))
    return Ldag

def mpow(R, p):
    w, V = np.linalg.eigh(R); return (V*np.clip(w, 1e-300, None)**p)@V.conj().T

def tracenorm(M): return np.abs(np.linalg.eigvalsh((M+M.conj().T)/2)).sum()

def run(name, n, g, beta=1.5, A=(0,), ts=(1, 10, 100, 1000, 10000)):
    d = 2**n; H = ham(n, g=g)
    rho = mpow(np.eye(d), 1)*0; w, V = np.linalg.eigh(H); rho = (V*np.exp(-beta*w))@V.conj().T; rho /= np.trace(rho)
    # rho_{-A} = tau_A (x) Tr_A rho, via the full Pauli twirl on A
    paulis1 = [I2, X, Y, Z]
    import itertools
    PA = []
    for combo in itertools.product(range(4), repeat=len(A)):
        S = np.eye(d, dtype=complex)
        for site, c in zip(A, combo): S = S@op(paulis1[c], site, n)
        PA.append((combo, S))
    rho_mA = sum(S@rho@S for _, S in PA)/4**len(A)
    # (a) twirl identity over non-trivial strings
    lhs = rho - rho_mA
    rhs = sum(S@(S@rho - rho@S) - (S@rho - rho@S)@S for combo, S in PA if any(combo))/2**(2*len(A)+1)
    print(f"[{name}] (a) twirl identity error {np.abs(lhs-rhs).max():.1e}")
    # (b) index bound
    print(f"[{name}] (b) min eig of 4^|A| rho_-A - rho = {np.linalg.eigvalsh(4**len(A)*rho_mA - rho).min():.2e} (>= 0 expected)")
    # generator
    jumps = [op(P, s, n) for s in A for P in (X, Y, Z)]
    Ldag = davies_heis(H, jumps, beta)
    S = np.kron(mpow(rho, 0.25), mpow(rho, 0.25).T)          # X -> rho^{1/4} X rho^{1/4}
    Si = np.kron(mpow(rho, -0.25), mpow(rho, -0.25).T)
    K = S@Ldag@Si
    print(f"[{name}] (c) KMS symmetry defect ||K - K^H|| = {np.abs(K-K.conj().T).max():.1e}")
    K = (K+K.conj().T)/2; k, Vk = np.linalg.eigh(K)
    zero = np.abs(k) < 1e-9
    Fdim = int(zero.sum())
    ker = Si@Vk[:, zero]                                      # basis of F_A = ker L^dag (Heisenberg operators)
    # N_A = 1_A (x) B(H_{A^c}) : orthonormal vec basis via the projector X -> (1/4^|A|) sum_S S X S
    EN = sum(np.kron(Sm, Sm.T) for _, Sm in PA)/4**len(A)     # tracial conditional expectation onto N_A (row-major vec)
    inN = np.abs(EN@ker - ker).max()
    print(f"[{name}] (d) dim F_A = {Fdim}, dim N_A = {4**(n-len(A))}, max ||E_N(x)-x|| over F_A basis = {inN:.1e}")
    # (e) exact recovery at t = infinity
    Pdag = Si@Vk[:, zero]@Vk[:, zero].conj().T@S             # KMS-orthogonal projection onto F_A
    P = Pdag.conj().T                                         # Schrodinger picture
    rec = (P@rho_mA.reshape(-1)).reshape(d, d)
    print(f"[{name}] (e) ||P_A[rho_-A] - rho||_1 = {tracenorm(rec-rho):.1e}")
    # (f) largest ad_H-invariant subspace of N_A
    adH = np.kron(H, np.eye(d)) - np.kron(np.eye(d), H.T)
    wN, VN = np.linalg.eigh((EN+EN.conj().T)/2); B = VN[:, wN > 0.5]
    while True:
        Pk = B@B.conj().T
        M = (np.eye(d*d) - Pk)@adH@B
        u, s, vh = np.linalg.svd(M, full_matrices=True)
        null = vh[np.sum(s > 1e-8):].conj().T
        Bn = B@null
        if Bn.shape[1] == B.shape[1]: break
        B = Bn
    PNs = B@B.conj().T
    print(f"[{name}] (f) dim N_A^sigma = {B.shape[1]}; max ||(1-P_Nsigma) x|| over F_A basis = {np.abs(ker - PNs@ker).max():.1e}")
    # (g) locality of P_A: action on an operator at the site farthest from A
    far = n//2
    Xfar = op(X, far, n)
    out = (Pdag@Xfar.reshape(-1)).reshape(d, d)
    print(f"[{name}] (g) ||P_A^dag(X_far) - X_far|| = {np.abs(out-Xfar).max():.2e},  ||P_A^dag(X_far) - Tr(rho X_far) 1|| = {np.abs(out-np.trace(rho@Xfar)*np.eye(d)).max():.2e}")
    # (h) finite-time recovery error
    gap = np.min(np.abs(k[~zero]))
    errs = []
    for t in ts:
        phi = np.where(zero, 1.0, (np.exp(t*k)-1)/(t*np.where(zero, 1, k)))
        Rdag = Si@(Vk*phi)@Vk.conj().T@S
        R = Rdag.conj().T
        errs.append(tracenorm((R@rho_mA.reshape(-1)).reshape(d, d) - rho))
    print(f"[{name}] (h) local gap above ker = {gap:.3e}; ||R_t[rho_-A]-rho||_1 at t={list(ts)}: {np.array(errs)}")

def ptrace_keep(R, keep, n):
    T = R.reshape([2]*(2*n)); idx = list(range(n)); out = list(range(n, 2*n))
    letters = 'abcdefghijklmnopqrstuvwxyz'
    a = [letters[i] for i in range(n)]; b = [letters[n+i] if i in keep else letters[i] for i in range(n)]
    expr = ''.join(a)+''.join(b)+'->'+''.join(a[i] for i in keep)+''.join(b[i] for i in keep)
    k = len(keep); return np.einsum(expr, T).reshape(2**k, 2**k)

def vn(R):
    w = np.linalg.eigvalsh((R+R.conj().T)/2); w = w[w > 1e-15]; return -(w*np.log(w)).sum()

def relent(R, S):
    def lg(M):
        w, V = np.linalg.eigh((M+M.conj().T)/2); return (V*np.log(w))@V.conj().T
    return np.trace(R@(lg(R)-lg(S))).real

def cmi_check(n=5, g=0.9, beta=1.5):
    d = 2**n; H = ham(n, g=g); w, V = np.linalg.eigh(H); rho = (V*np.exp(-beta*w))@V.conj().T; rho /= np.trace(rho)
    A, B, C = [0], [1, n-1], [k for k in range(2, n-1)]
    S = lambda keep: vn(ptrace_keep(rho, sorted(keep), n))
    cmi = S(A+B) + S(B+C) - S(B) - S(A+B+C)
    rho_mA = sum(op(P, 0, n)@rho@op(P, 0, n) for P in (I2, X, Y, Z))/4
    delta = relent(rho, rho_mA) - relent(ptrace_keep(rho, sorted(A+B), n), ptrace_keep(rho_mA, sorted(A+B), n))
    print(f"[TFIM n={n}] (i) I(A:C|B) = {cmi:.6e},  delta_AB(rho, rho_-A) = {delta:.6e},  diff = {abs(cmi-delta):.1e}")
    # (j) reflection i -> -i mod n fixes site 0
    perm = [(-i) % n for i in range(n)]
    R = np.zeros((d, d))
    for s in range(d):
        bits = [(s >> (n-1-i)) & 1 for i in range(n)]
        nb = [bits[perm[i]] for i in range(n)]
        R[int(''.join(map(str, nb)), 2), s] = 1
    jumps = [op(P, 0, n) for P in (X, Y, Z)]
    Ldag = davies_heis(H, jumps, beta)
    print(f"[TFIM n={n}] (j) ||[H,R]|| = {np.abs(H@R-R@H).max():.1e}, ||L_A^dag(R)|| = {np.abs(Ldag@R.reshape(-1).astype(complex)).max():.1e}")

if __name__ == "__main__":
    cmi_check()
    run("classical Ising ring, n=5", n=5, g=0.0)
    run("transverse-field Ising ring, n=5", n=5, g=0.9)
