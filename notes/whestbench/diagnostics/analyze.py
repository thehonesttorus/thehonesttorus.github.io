# Closure diagnostics from pass2_stats.npz (network 0). Oracle inputs: mu from 1e9 truth means, Sigma from 1.5e6 samples.
import numpy as np, time, sys
from scipy.special import ndtr
t0=time.time()
W=np.load('data/W0.npy').astype(np.float64); Y=np.load('data/Y0.npy').astype(np.float64); n=1024; L=16
S=np.load('pass2_stats.npz'); Hm=S['Hm']; HS=S['HS']; Sz=S['Sz']; N2=int(S['N2'])
SQ2PI=np.sqrt(2*np.pi)
def pdf(t): return np.exp(-0.5*t*t)/SQ2PI
def He(k,x):
    if k==0: return np.ones_like(x)
    if k==1: return x
    a,b=np.ones_like(x),x
    for m in range(1,k): a,b=b,x*b-m*a
    return b
def dcoef(k,t):   # Hermite coefficients of (t+zeta)_+ : d_1=Phi(t), d_k=phi(t)He_{k-2}(-t)
    return ndtr(t) if k==1 else pdf(t)*He(k-2,-t)
def relu_moments(mu,s):   # E[z^k 1{z>0}], k=0..4
    t=mu/s; m=[ndtr(t)]; g0=pdf(t)/s
    m.append(mu*m[0]+s*s*g0)
    for k in range(2,5): m.append(mu*m[k-1]+s*s*(k-1)*m[k-2])
    return m
def cums(M):
    a,b,c,d=M; v=b-a*a; k3=c-3*a*b+2*a**3; k4=d-4*a*c-3*b*b+12*a*a*b-6*a**4
    return a,v,k3,k4
def edge(mu,s,k3,k4,order):
    t=mu/s; g3=k3/s**3; g4=k4/s**4
    corr=(g3/6)*(-t)
    if order>=4: corr=corr+(g4/24)*(t*t-1)+(g3*g3/72)*(t**4-6*t*t+3)
    return s*pdf(t)*corr
def mse(a,b): return float(np.mean((a-b)**2))
# per-layer oracle inputs
mu=[None]*L; Szz=[None]*L; sig=[None]*L; tt=[None]*L
for l in range(L):
    mu[l]=(Y[l-1]@W[l]) if l>0 else np.zeros(n)
    Szz[l]=W[l].T@HS[l]@W[l]; sig[l]=np.sqrt(np.diag(Szz[l])); tt[l]=mu[l]/sig[l]
# source models (third cumulant of h_{l+1}=phi(z_l) at leading order): M=D1 Sz (off-diag), p=phi(t)/s, q=2Psi(1-Phi), k3h
def source(l):
    m=relu_moments(mu[l],sig[l]); d1=ndtr(tt[l]); Psi=m[1]
    M=d1[:,None]*Szz[l]; np.fill_diagonal(M,0.0)
    p=pdf(tt[l])/sig[l]; q=2*Psi*(1-d1)
    k3h=m[3]-3*m[2]*m[1]+2*m[1]**3
    return M,p,q,k3h,d1
src=[source(l) for l in range(L)]
def k3_model(l,K):
    """modelled third cumulant of z_l from sources s=l-1,...,l-K (back-propagated directions)."""
    out=np.zeros(n); U=W[l].copy()      # columns: directions in h_l space
    for s in range(l-1,max(l-1-K,-1),-1):
        # U currently lives in h_{s+1} space
        M,p,q,k3h,_=src[s]
        Nt=M.T@U                          # (n x n): (M^T u_j)_b
        out+=(3*p[:,None]*U*Nt*Nt+3*q[:,None]*U*U*Nt+k3h[:,None]*U**3).sum(0)
        if s>0: U=W[s]@(src[s-1][4][:,None]*U)   # back-propagate through layer s: u <- W[s] D1^{(s-1)} u ... 
    return out
# NB: h_{s+1}=phi(z_s), z_s=h_s@W[s]; linear response dz_s = W[s]^T dh_s, dh_{s+1}=D1^{(s)} dz_s -> direction map u_{h_s}=W[s] D1^{(s)} u_{h_{s+1}}
def k3_model_fixed(l,K):
    out=np.zeros(n); U=W[l].copy()
    for s in range(l-1,max(l-1-K,-1),-1):
        M,p,q,k3h,d1=src[s]
        Nt=M.T@U
        out+=(3*p[:,None]*U*Nt*Nt+3*q[:,None]*U*U*Nt+k3h[:,None]*U**3).sum(0)
        U=W[s]@(d1[:,None]*U)
    return out
print("layer |  G       E3(smp)  E4(smp) | E3 with modelled k3: K=1    K=2    K=3    K=5    K=all | R2(k3 model vs sample) K=1 K=2 K=3 K=5 all | skew rms  kurt rms")
for l in range(1,L):
    truth=Y[l]; a,v,k3,k4=cums(Sz[l]); s=sig[l]; m=mu[l]
    G=m*ndtr(m/s)+s*pdf(m/s)
    E3=G+edge(m,s,k3,k4,3); E4=G+edge(m,s,k3,k4,4)
    row=[mse(G,truth),mse(E3,truth),mse(E4,truth)]
    r2=[]; e3m=[]
    for K in [1,2,3,5,l]:
        km=k3_model_fixed(l,K)
        r2.append(1-np.mean((km-k3)**2)/np.mean(k3**2))
        e3m.append(mse(G+edge(m,s,km,k4,3),truth))
    g3=k3/s**3; g4=k4/s**4
    print(f"{l+1:5d} | {row[0]:.1e} {row[1]:.1e} {row[2]:.1e} | {e3m[0]:.1e} {e3m[1]:.1e} {e3m[2]:.1e} {e3m[3]:.1e} {e3m[4]:.1e} | {r2[0]:.3f} {r2[1]:.3f} {r2[2]:.3f} {r2[3]:.3f} {r2[4]:.3f} | {np.sqrt(np.mean(g3**2)):.3f} {np.sqrt(np.mean(g4**2)):.3f}", flush=True)
print(f"[{time.time()-t0:.0f}s] covariance closure check: Gaussian closure cov of h_{{l+1}} from oracle (mu,Sz) vs sample cov")
print("layer | diag rel err mean, rms | next-layer sigma^2 rel err mean, rms | offdiag-only contribution to next sigma^2: rel err mean, rms | max|rho|")
for l in range(1,L-1):
    s=sig[l]; t=tt[l]; Sz_=Szz[l]; rho=Sz_/np.outer(s,s); np.fill_diagonal(rho,0)
    C=np.zeros((n,n)); rk=np.ones((n,n))
    for k in range(1,13):
        rk=rk*rho; dk=dcoef(k,t); C+=rk*np.outer(dk,dk)/np.math.factorial(k) if hasattr(np,'math') else rk*np.outer(dk,dk)/float(__import__('math').factorial(k))
    C*=np.outer(s,s)
    mm=relu_moments(mu[l],s); var_exact=mm[2]-mm[1]**2
    np.fill_diagonal(C,var_exact)
    Strue=HS[l+1]
    dre=(np.diag(C)-np.diag(Strue))/np.diag(Strue)
    Wn=W[l+1]; v_model=np.einsum('ij,ik,kj->j',Wn,C,Wn); v_true=np.einsum('ij,ik,kj->j',Wn,Strue,Wn)
    re=(v_model-v_true)/v_true
    Coff=C.copy(); np.fill_diagonal(Coff,0); Soff=Strue.copy(); np.fill_diagonal(Soff,0)
    reo=(np.einsum('ij,ik,kj->j',Wn,Coff,Wn)-np.einsum('ij,ik,kj->j',Wn,Soff,Wn))/v_true
    print(f"{l+1:5d} | {dre.mean():+.2e} {np.sqrt(np.mean(dre**2)):.2e} | {re.mean():+.2e} {np.sqrt(np.mean(re**2)):.2e} | {reo.mean():+.2e} {np.sqrt(np.mean(reo**2)):.2e} | {np.abs(rho).max():.3f}", flush=True)
print(f"total {time.time()-t0:.0f}s")
