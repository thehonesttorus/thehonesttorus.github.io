# Section 5 claim: a covariance-matched trace-cubic source at a CENTERED reference is exactly invisible
# to every degree-one homogeneous observable (so to every bias-free ReLU suffix mean), but it is visible
# to intermediate second moments.  A Gaussian closure that keeps the second-moment change and drops the
# compensating structure can therefore manufacture a mean response out of a provably null sector.
#
# Truth side: the first-order response is  delta_R M = E_gamma[ s(Y) f(Sigma^{1/2} Y) ],
#   s(y) = (1/6) R(w,I):H_3(y) = (1/2)(w.y)(||y||^2 - (n+2)),
# which is exactly zero by radial cancellation: E[R^2(R^2-(n+2))] = 0 for R ~ chi_n.
# Closure side: inject the same tensor into our propagation and see what it returns.
import numpy as np, time
from scipy.special import ndtr
import twin2 as T
n=1024; L=16; pdf=T.pdf
seeds,Yall=T.load_nets(); W=[w.astype(np.float64) for w in T.regen(seeds[0])]
HS=np.load('../pass2_stats.npz')['HS']
def psqrt(S):
    ev,V=np.linalg.eigh(S); ev=np.clip(ev,0,None); return (V*np.sqrt(ev))@V.T
LAYER=8
Sig=HS[LAYER].copy(); Sig=Sig+1e-3*np.mean(np.diag(Sig))*np.eye(n)
A=psqrt(Sig); Ainv=np.linalg.pinv(A)
rs=np.random.default_rng(3); v=rs.standard_normal(n); v/=np.linalg.norm(v)
w=Ainv@v
# ---------- truth: the exact first-order response of the suffix mean ----------
def suffix(H,l):
    for k in range(l,L): H=np.maximum(H@W[k],0.0)
    return H
NS=400_000; CH=20_000; rng=np.random.default_rng(2718)
acc=np.zeros(n); acc0=np.zeros(n); accs=0.0; N=0; a2=np.zeros(n)
for _ in range(NS//CH):
    y=rng.standard_normal((CH,n))
    s=0.5*(y@w)*((y*y).sum(1)-(n+2))
    H=suffix(y@A, LAYER)
    acc+=(s[:,None]*H).sum(0); acc0+=H.sum(0); a2+=((s[:,None]*H)**2).sum(0); N+=CH
resp=acc/N; base=acc0/N
se=np.sqrt(np.maximum(a2/N-resp**2,0)/N)
print(f"TRUTH, centered reference at layer {LAYER}, suffix -> layer 16, {NS} samples")
print(f"  ||delta_R M||          = {np.linalg.norm(resp):.4e}")
print(f"  Monte-Carlo noise floor = {np.linalg.norm(se):.4e}   (identity predicts exactly 0)")
print(f"  ||M||                   = {np.linalg.norm(base):.4e}   ratio {np.linalg.norm(resp)/np.linalg.norm(base):.2e}")
# ---------- closure: inject the same tensor and propagate with our equations ----------
def closure(inject_eps):
    m=np.zeros(n); Sg=Sig.copy(); outs=[]
    for l in range(LAYER,L):
        Wl=W[l]; mu=m@Wl; Sz=(Wl.T@Sg)@Wl; sz2=np.diag(Sz).copy(); sz=np.sqrt(sz2); t=mu/sz
        ph=pdf(t); Ph=ndtr(t); Psi=mu*Ph+sz*ph; m2=(mu*mu+sz2)*Ph+mu*sz*ph
        s3=sz2*sz; s4=sz2*sz2
        k3=np.zeros(n); K21=np.zeros((n,n))
        if l==LAYER and inject_eps!=0.0:
            vt=Wl.T@v                                   # R(v,Sigma) transports to R(W^T v, W^T Sigma W)
            k3=inject_eps*3.0*vt*sz2                    # kappa_3(z)_{aaa}
            K21=inject_eps*(2.0*vt[:,None]*Sz+vt[None,:]*sz2[:,None])   # kappa_3(z)_{aab}
            np.fill_diagonal(K21,0.0)
        t2=t*t; he4=(t2-6)*t2+3
        dE1=-k3*t*ph/(6*sz2)+(k3*k3)*he4*ph/(72*s4*sz)
        dE2=k3*ph/(3*sz)+(k3*k3)*(3*t-t*t2)*ph/(36*s4)
        m_next=Psi+dE1
        so=np.outer(sz,sz); rho=Sz/so; ds=T.hermite_d(t,ph,Ph,8); rk=rho.copy(); C=np.zeros((n,n)); fact=1.0
        for kk in range(1,9):
            fact*=kk; C+=rk*np.outer(ds[kk-1],ds[kk-1])/fact
            if kk<8: rk*=rho
        C*=so
        F2=ph/sz
        C=C+0.5*(K21*np.outer(F2,Ph)+K21.T*np.outer(Ph,F2))
        C[np.diag_indices(n)]=(m2+dE2)-m_next*m_next
        m=m_next; Sg=C; outs.append(m)
    return outs[-1]
print("\nCLOSURE response to the same injected tensor (truth says 0):")
b0=closure(0.0)
for eps in (1e-3,3e-3,1e-2):
    d=closure(eps)-b0
    print(f"  eps={eps:.0e}: ||closure response||/eps = {np.linalg.norm(d)/eps:.4e}"
          f"   relative to ||M|| = {np.linalg.norm(d)/eps/np.linalg.norm(base):.2e}")

# ---- sizing: is the spurious response comparable to the closure's response to a legitimate cubic
# of the same norm?  Inject a norm-matched generic (non-trace) k3/K21 pair and compare.
def closure_generic(eps, seed):
    m=np.zeros(n); Sg=Sig.copy(); outs=[]
    rr=np.random.default_rng(seed)
    for l in range(LAYER,L):
        Wl=W[l]; mu=m@Wl; Sz=(Wl.T@Sg)@Wl; sz2=np.diag(Sz).copy(); sz=np.sqrt(sz2); t=mu/sz
        ph=pdf(t); Ph=ndtr(t); Psi=mu*Ph+sz*ph; m2=(mu*mu+sz2)*Ph+mu*sz*ph
        s3=sz2*sz; s4=sz2*sz2
        k3=np.zeros(n); K21=np.zeros((n,n))
        if l==LAYER and eps!=0.0:
            vt=Wl.T@v
            k3r=rr.standard_normal(n); k3r*=np.linalg.norm(3.0*vt*sz2)/np.linalg.norm(k3r)
            Kr=rr.standard_normal((n,n)); np.fill_diagonal(Kr,0.0)
            Kt=2.0*vt[:,None]*Sz+vt[None,:]*sz2[:,None]; np.fill_diagonal(Kt,0.0)
            Kr*=np.linalg.norm(Kt)/np.linalg.norm(Kr)
            k3=eps*k3r; K21=eps*Kr
        t2=t*t; he4=(t2-6)*t2+3
        dE1=-k3*t*ph/(6*sz2)+(k3*k3)*he4*ph/(72*s4*sz)
        dE2=k3*ph/(3*sz)+(k3*k3)*(3*t-t*t2)*ph/(36*s4)
        m_next=Psi+dE1
        so=np.outer(sz,sz); rho=Sz/so; ds=T.hermite_d(t,ph,Ph,8); rk=rho.copy(); C=np.zeros((n,n)); fact=1.0
        for kk in range(1,9):
            fact*=kk; C+=rk*np.outer(ds[kk-1],ds[kk-1])/fact
            if kk<8: rk*=rho
        C*=so; F2=ph/sz
        C=C+0.5*(K21*np.outer(F2,Ph)+K21.T*np.outer(Ph,F2))
        C[np.diag_indices(n)]=(m2+dE2)-m_next*m_next
        m=m_next; Sg=C; outs.append(m)
    return outs[-1]
print("\nSizing: closure response per unit injection norm, trace sector vs norm-matched generic")
b0=closure(0.0); eps=1e-3
dt=np.linalg.norm(closure(eps)-b0)/eps
gs=[np.linalg.norm(closure_generic(eps,sd)-b0)/eps for sd in (11,12,13)]
print(f"  trace sector R(v,Sigma)      : {dt:.4e}")
print(f"  norm-matched generic (3 draws): {gs[0]:.4e} {gs[1]:.4e} {gs[2]:.4e}   mean {np.mean(gs):.4e}")
print(f"  ratio trace/generic = {dt/np.mean(gs):.3f}")
