# The Mac audit's freezing result is about a MOVING selection: P(s) A(s) P(s), where d/ds mixes the
# change of the state with the change of the chosen resolution.  Our estimator has exactly one such
# object -- T_ACT = 2.75 gates which neurons get the exact two-point treatment, by |t| < 2.75, and t
# moves with the state.  Two questions: how much does the gate actually exclude, and does it matter?
import numpy as np
from scipy.special import ndtr
import ablate as AB, twin2 as T, cov_rho as CR
n=1024; L=16; INV=T.INV; pdf=T.pdf
A_S,G_S,C_S=AB.A_S,AB.G_S,AB.C_S
def gate_stats(k,N=2400,seed=0):
    W=AB.getW(k)
    rng=np.random.default_rng(seed); X=rng.standard_normal((N,n)).astype(np.float32); Zs=[]; h=X
    for l in range(L): z=h@W[l]; Zs.append(z); h=np.maximum(z,0)
    sources={}; out=[]
    W0=W[0]; G=W0.T@W0; s2=np.diag(G); s=np.sqrt(s2); m=s*INV; so0=np.outer(s,s)
    rho=np.clip(G/so0,-1,1); C=so0*((np.sqrt(1-rho*rho)+rho*(np.pi-np.arccos(rho)))/(2*np.pi))-np.outer(m,m)
    C[np.diag_indices(n)]=s2*(0.5-1/(2*np.pi)); C=C.astype(np.float32)
    M=(G*0.5).astype(np.float32); M[np.diag_indices(n)]=0
    sources[0]=(M,(INV/s).astype(np.float32),m.astype(np.float32),
                (s2*s*(2*INV-1.5*INV+2*INV**3)).astype(np.float32),np.full(n,0.5,np.float32))
    Sig=C; m=m.astype(np.float32); prev=None
    for l in range(1,L):
        Wl=W[l]; mu=(m@Wl).astype(np.float64); Sz=((Wl.T@Sig)@Wl).astype(np.float64)
        sz2=np.diag(Sz).copy(); sz=np.sqrt(np.maximum(sz2,1e-30)); t=mu/sz
        act=np.abs(t)<2.75
        out.append((l,int(act.sum()),float(np.abs(t).max()),
                    -1 if prev is None else int((act!=prev).sum())))
        prev=act
        ph=pdf(t); Ph=ndtr(t); Psi=mu*Ph+sz*ph; m2=(mu*mu+sz2)*Ph+mu*sz*ph; s3=sz2*sz; s4=sz2*sz2
        Ms,ps,qs,k3hs,d1s=sources[l-1]; Nt=Ms.T@Wl; UN=Wl*Nt
        k3=((ps[:,None]*np.float32(3))*UN*Nt+(qs[:,None]*np.float32(3))*UN*Wl+k3hs[:,None]*(Wl*Wl*Wl)).sum(0)
        X1=(UN*2)*ps[:,None]+(Wl*Wl)*qs[:,None]
        X2=(Nt*Nt)*ps[:,None]+(UN*2)*qs[:,None]+(Wl*Wl)*k3hs[:,None]
        K21=(X1.T@Nt+X2.T@Wl).astype(np.float64); np.fill_diagonal(K21,0.0)
        k3=k3.astype(np.float64)+A_S[l-1]*t*s3; g4=G_S[l-1]
        z=Zs[l]; a=z.mean(0); zc=z-a; zc2=zc*zc; k3s=(zc2*zc).mean(0).astype(np.float64)
        rr=k3s-k3; nz=6*np.mean(s3*s3)/N; sg=max(np.mean(rr**2)-nz,0); k3=k3+(sg/(sg+nz))*rr
        muh=mu/np.linalg.norm(mu); Ap=zc@muh.astype(np.float32)
        u=0.75*(C_S[l-1]*s3)+0.25*((zc2.T@Ap).astype(np.float64)/N-K21@muh)
        K21=K21+np.outer(u,muh); np.fill_diagonal(K21,0.0)
        k4v=g4*s4; t2=t*t; he2=t2-1; he4=(t2-6)*t2+3
        dE1=-k3*t*ph/(6*sz2)+(k3*k3)*he4*ph/(72*s4*sz)+k4v*he2*ph/(24*s3)
        dE2=k3*ph/(3*sz)+(k3*k3)*(3*t-t*t2)*ph/(36*s4)-k4v*t*ph/(12*sz2)
        m_next=Psi+dE1
        so=np.outer(sz,sz); rhoc=Sz/so; ds=T.hermite_d(t,ph,Ph,8); rk=rhoc.copy()
        Cn=np.zeros((n,n),np.float32); fact=1.0
        for kk in range(1,9):
            fact*=kk; Cn=Cn+(rk*np.outer(ds[kk-1],ds[kk-1])/fact).astype(np.float32)
            if kk<8: rk*=rhoc
        Cn=Cn*so.astype(np.float32); F2=ph/sz
        Cn=Cn+(0.5*(K21*np.outer(F2,Ph)+K21.T*np.outer(Ph,F2))).astype(np.float32)
        v4=sz*ph; tp=t*ph
        Cn=Cn+(0.0020*np.outer(v4,v4)).astype(np.float32)
        Cn=Cn+(0.22*g4*Sz*(np.outer(tp,Ph)+np.outer(Ph,tp))).astype(np.float32)
        Cn[np.diag_indices(n)]=((m2+dE2)-m_next*m_next).astype(np.float32)
        Msrc=Sz*Ph[:,None]; np.fill_diagonal(Msrc,0.0); m3=mu*m2+2*sz2*Psi
        sources[l]=(Msrc.astype(np.float32),(ph/sz).astype(np.float32),(2*Psi*(1-Ph)).astype(np.float32),
                    (m3-3*m2*Psi+2*Psi**3).astype(np.float32),Ph.astype(np.float32))
        Sig=Cn.astype(np.float32); m=m_next.astype(np.float32)
    return out
for k in (0,1):
    print(f"network {k}:  |t| < 2.75 gate, out of {n} neurons")
    print(f"  {'layer':>5} {'active':>7} {'excluded':>9} {'max|t|':>8} {'membership churn vs prev':>25}")
    for l,a,mx,ch in gate_stats(k):
        print(f"  {l:5d} {a:7d} {n-a:9d} {mx:8.3f} {('-' if ch<0 else str(ch)):>25}")
    print()
