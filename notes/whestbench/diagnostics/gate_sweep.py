# The gate is a real moving selection: 40% excluded by depth, ~500 neurons changing membership per
# layer.  It does not corrupt the VALUE (it is just which rows of K21 get computed), but it is worth
# asking what it costs in accuracy and what widening it costs in budget.
# Cost model (section 32): the two-point contraction is 4.0 * (n_act/n) n^3 per layer; the rest of the
# estimator is 198 - 15*3.2 = 150 n^3.  Budget 2048 n^3, floor at 10% = 204.8 n^3.
import numpy as np
from scipy.special import ndtr
import ablate as AB, twin2 as T
n=1024; L=16; INV=T.INV; pdf=T.pdf
A_S,G_S,C_S=AB.A_S,AB.G_S,AB.C_S
def predict(k,tact=np.inf,N=2400,seed=0):
    W=AB.getW(k); Y=AB.Yall[k]
    rng=np.random.default_rng(seed); X=rng.standard_normal((N,n)).astype(np.float32); Zs=[]; h=X
    for l in range(L): z=h@W[l]; Zs.append(z); h=np.maximum(z,0)
    sources={}; acts=[]
    W0=W[0]; G=W0.T@W0; s2=np.diag(G); s=np.sqrt(s2); m=s*INV; so0=np.outer(s,s)
    rho=np.clip(G/so0,-1,1); C=so0*((np.sqrt(1-rho*rho)+rho*(np.pi-np.arccos(rho)))/(2*np.pi))-np.outer(m,m)
    C[np.diag_indices(n)]=s2*(0.5-1/(2*np.pi)); C=C.astype(np.float32)
    M=(G*0.5).astype(np.float32); M[np.diag_indices(n)]=0
    sources[0]=(M,(INV/s).astype(np.float32),m.astype(np.float32),
                (s2*s*(2*INV-1.5*INV+2*INV**3)).astype(np.float32),np.full(n,0.5,np.float32))
    Sig=C; m=m.astype(np.float32)
    for l in range(1,L):
        Wl=W[l]; mu=(m@Wl).astype(np.float64); Sz=((Wl.T@Sig)@Wl).astype(np.float64)
        sz2=np.diag(Sz).copy(); sz=np.sqrt(np.maximum(sz2,1e-30)); t=mu/sz
        ph=pdf(t); Ph=ndtr(t); Psi=mu*Ph+sz*ph; m2=(mu*mu+sz2)*Ph+mu*sz*ph; s3=sz2*sz; s4=sz2*sz2
        Ms,ps,qs,k3hs,d1s=sources[l-1]; Nt=Ms.T@Wl; UN=Wl*Nt
        k3=((ps[:,None]*np.float32(3))*UN*Nt+(qs[:,None]*np.float32(3))*UN*Wl+k3hs[:,None]*(Wl*Wl*Wl)).sum(0)
        X1=(UN*2)*ps[:,None]+(Wl*Wl)*qs[:,None]
        X2=(Nt*Nt)*ps[:,None]+(UN*2)*qs[:,None]+(Wl*Wl)*k3hs[:,None]
        act=np.abs(t)<tact if np.isfinite(tact) else np.ones(n,bool)
        acts.append(act.mean())
        K21=np.zeros((n,n))
        Ka=(X1[:,act].T@Nt+X2[:,act].T@Wl).astype(np.float64)
        K21[act,:]=Ka; np.fill_diagonal(K21,0.0)
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
    return float(np.mean((m.astype(np.float64)-Y[-1])**2)), float(np.mean(acts))
NETS=list(range(16))
print("T_ACT sweep on the 16 local MLPs.  cost = 150 + 15*4.0*frac_active (n^3); budget 2048;")
print("the 10% floor is 204.8 n^3, so anything above that is charged linearly.\n")
print(f"{'T_ACT':>7} {'frac active':>12} {'cost n^3':>9} {'C/B':>7} {'mult':>6} {'raw MSE':>11} {'adjusted':>11}")
for ta in (2.75, 3.25, 3.75, 4.5, np.inf):
    res=[predict(k,ta) for k in NETS]
    mse=float(np.mean([r[0] for r in res])); fr=float(np.mean([r[1] for r in res]))
    cost=150.0+15*4.0*fr; cb=cost/2048.0; mult=max(0.1,cb)
    print(f"{ta:7.2f} {fr:12.4f} {cost:9.1f} {cb:7.2%} {mult:6.3f} {mse:11.4e} {mse*mult:11.4e}", flush=True)
