# Own numerical check (not from any source): Lemma III.3 identities and the bound j(s) <= (1+s) q^2 / s^2
# of arXiv:2609.38007 for a general subalgebra N = (+) M_k (x) 1_m with trace-preserving E_N. Run: python3 this_file.py
import numpy as np
exec(open(__import__('os').path.join(__import__('os').path.dirname(__file__), 'check_2609_subalgebra_bound.py')).read().split("for blocks in")[0])   # reuse helpers (rand_state, make_EN, fpow, fexpit, relent)
rng2 = np.random.default_rng(7)
def vec(X): return X.reshape(-1)
def L(A): n=A.shape[0]; return np.kron(A, np.eye(n))          # X -> A X  (row-major vec)
def R(A): n=A.shape[0]; return np.kron(np.eye(n), A.T)        # X -> X A
worst = {}
for blocks in [[(2,2),(1,3)], [(3,2)], [(1,1),(2,1),(1,2)]]:
    n, EN = make_EN(blocks)
    rho, sig = rand_state(n, 2.5), rand_state(n, 2.5)
    rN, sN = EN(rho), EN(sig)
    # basis of N as operators: images of E_N on matrix units, orthonormalised in HS
    imgs = []
    for i in range(n):
        for j in range(n):
            E = np.zeros((n,n),complex); E[i,j]=1; imgs.append(vec(EN(E)))
    Q, Rr = np.linalg.qr(np.array(imgs).T); rank = np.sum(np.abs(np.diag(Rr))>1e-10)
    U, S, _ = np.linalg.svd(np.array(imgs).T); BN = U[:, :int(np.sum(S>1e-10))]   # ONB of N in HS space
    W = np.column_stack([vec(b.reshape(n,n) @ fpow(rN,-0.5) @ fpow(rho,0.5)) for b in BN.T])
    Delta = L(sig) @ R(np.linalg.inv(rho)); DeltaN = BN.conj().T @ (L(sN) @ R(np.linalg.inv(rN))) @ BN
    err_iso = np.linalg.norm(W.conj().T@W - np.eye(W.shape[1]))
    err_comp = np.linalg.norm(W.conj().T@Delta@W - DeltaN)
    # j(s) and the bound (1+s) q^2 / s^2
    ts = np.linspace(-12,12,4801); ts=ts[np.abs(ts)>1e-9]; dt=ts[1]-ts[0]
    eta = np.array([np.linalg.norm((lambda u: u-EN(u))(fexpit(sig,t)@fexpit(rho,-t)),2) for t in ts])
    q = 0.5*np.sum(eta/np.abs(np.sinh(np.pi*ts)))*dt
    r12 = vec(fpow(rho,0.5)); rN12 = BN.conj().T @ vec(fpow(rN,0.5))
    def j(s):
        return (r12.conj() @ np.linalg.solve(Delta + s*np.eye(n*n), r12)).real - (rN12.conj() @ np.linalg.solve(DeltaN + s*np.eye(len(rN12)), rN12)).real
    ratios = [j(s)/((1+s)*q*q/s**2) for s in np.logspace(-3,3,25)]
    # delta via resolvent integral vs direct
    ss = np.logspace(-6,6,4000); js = np.array([j(s) for s in ss]); integ = np.trapezoid(js, ss)
    delta = relent(rho,sig) - relent(rN,sN)
    print(blocks, f"iso_err={err_iso:.1e} WdagDeltaW-DeltaN={err_comp:.1e} max j/bound={max(ratios):.3f} min j={min(js):.1e} delta={delta:.5f} int j={integ:.5f}")
