# One more round with S3 in place.  Candidates, all elementwise and all free:
#   U1 = Sz . (s phi) (x) (s phi)          the A4 term at next order in rho
#   U2 = Sz . Sz . [ (t phi) (x) Phi + T ] the S3 direction at rho^2
#   U3 = K21 . Sz . [ (t phi/s^2) (x) (phi/s) + T ]   the (a,a,b) term at next order in rho
#   U4 = Sz . [ (phi) (x) (phi) ]          the (a,a,b,b) pattern with a Phi-free weight
import numpy as np
from scipy.special import ndtr
import ablate as AB, twin2 as T
n=1024; L=16; INV=T.INV; pdf=T.pdf
seeds,Yall=AB.seeds,AB.Yall
A_S,G_S,C_S=AB.A_S,AB.G_S,AB.C_S
def predict(k,U=None,N=2400,seed=0,A4=0.0020,S3=0.22):
    U=U or {}
    W=AB.getW(k); Y=Yall[k]
    rng=np.random.default_rng(seed); X=rng.standard_normal((N,n)).astype(np.float32); Zs=[]; h=X
    for l in range(L): z=h@W[l]; Zs.append(z); h=np.maximum(z,0)
    sources={}
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
        Uw=Wl; Ms,ps,qs,k3hs,d1s=sources[l-1]; Nt=Ms.T@Uw; UN=Uw*Nt
        k3=((ps[:,None]*np.float32(3))*UN*Nt+(qs[:,None]*np.float32(3))*UN*Uw+k3hs[:,None]*(Uw*Uw*Uw)).sum(0)
        X1=(UN*2)*ps[:,None]+(Uw*Uw)*qs[:,None]
        X2=(Nt*Nt)*ps[:,None]+(UN*2)*qs[:,None]+(Uw*Uw)*k3hs[:,None]
        K21=(X1.T@Nt+X2.T@Uw).astype(np.float64); np.fill_diagonal(K21,0.0)
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
        Cn=Cn+(A4*np.outer(v4,v4)).astype(np.float32)
        Cn=Cn+(S3*g4*Sz*(np.outer(tp,Ph)+np.outer(Ph,tp))).astype(np.float32)
        if U.get('U1'): Cn=Cn+(U['U1']*Sz*np.outer(v4,v4)).astype(np.float32)
        if U.get('U2'): Cn=Cn+(U['U2']*g4*Sz*Sz*(np.outer(tp,Ph)+np.outer(Ph,tp))).astype(np.float32)
        if U.get('U3'):
            w1=t*ph/sz2; w2=ph/sz
            Cn=Cn+(U['U3']*Sz*(K21*np.outer(w1,w2)+K21.T*np.outer(w2,w1))).astype(np.float32)
        if U.get('U4'): Cn=Cn+(U['U4']*Sz*np.outer(ph,ph)).astype(np.float32)
        Cn[np.diag_indices(n)]=((m2+dE2)-m_next*m_next).astype(np.float32)
        Msrc=Sz*Ph[:,None]; np.fill_diagonal(Msrc,0.0); m3=mu*m2+2*sz2*Psi
        sources[l]=(Msrc.astype(np.float32),(ph/sz).astype(np.float32),(2*Psi*(1-Ph)).astype(np.float32),
                    (m3-3*m2*Psi+2*Psi**3).astype(np.float32),Ph.astype(np.float32))
        Sig=Cn.astype(np.float32); m=m_next.astype(np.float32)
    return float(np.mean((m.astype(np.float64)-Y[-1])**2))
if __name__=="__main__":
    FIT=[0,2,4,6]; HOLD=list(range(8,16))
    def sc(nets,U=None): return float(np.mean([predict(k,U) for k in nets]))
    b=sc(FIT); bh=sc(HOLD)
    print(f"shipped (A4 + S3):  fit {b:.4e}   held {bh:.4e}\n")
    print(f"{'spec':18s}  {'fit':>11} {'fit ratio':>10} {'held':>11} {'held ratio':>11}")
    for name,gs in (('U1',(-0.03,-0.01,0.01,0.03)),('U2',(-3.0,-1.0,1.0,3.0)),
                    ('U3',(-1.0,-0.3,0.3,1.0)),('U4',(-0.02,-0.006,0.006,0.02))):
        for v in gs:
            f=sc(FIT,{name:v})
            h=sc(HOLD,{name:v}) if f<b else float('nan')
            print(f"{name+'='+format(v,'.4g'):18s}  {f:11.4e} {f/b:10.4f} {h:11.4e} {h/bh:11.4f}", flush=True)
