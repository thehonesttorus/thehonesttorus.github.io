# The propagated pre-activation sigma is systematically too SMALL, monotonically in depth:
#   sigma_true/sigma_model - 1  =  6.0e-5, 1.0e-4, ... , 1.07e-3   (layers 1..15, nets 0 and 1, N=1e6)
# The Monte-Carlo error on that per-layer mean is 2.2e-5, so this is a real bias, not noise.  It comes
# from the OFF-DIAGONAL covariance closure: sigma^2(z_b) = sum_ac W_ab W_cb Sigma_ac is dominated by
# the n^2 off-diagonal terms, whose non-Gaussian correction we carry only to first order (the K21
# channel -- which is exactly the one section 32 found load-bearing).
# A per-layer multiplicative calibration of Sz costs n^2 flops, i.e. nothing.  Test it.
import numpy as np, sys
from scipy.special import ndtr
import ablate as AB, twin2 as T
n=1024; L=16; INV=T.INV; pdf=T.pdf
seeds,Yall=AB.seeds,AB.Yall
A_S,G_S,C_S=AB.A_S,AB.G_S,AB.C_S
E=np.array([5.97e-05,9.97e-05,8.89e-05,1.03e-04,1.74e-04,1.31e-04,2.06e-04,1.24e-04,2.29e-04,
            3.90e-04,5.57e-04,6.22e-04,7.07e-04,9.78e-04,1.07e-03])
V_CUM=2*E
V_INC=2*np.diff(np.concatenate([[0.0],E]))

def predict(k,V,N=2400,seed=0,al=1.0):
    W=AB.getW(k); Y=Yall[k]
    rng=np.random.default_rng(seed); X=rng.standard_normal((N,n)).astype(np.float32); Zs=[]; h=X
    for l in range(L): z=h@W[l]; Zs.append(z); h=np.maximum(z,0)
    sources={}
    W0=W[0]; G=W0.T@W0; s2=np.diag(G); s=np.sqrt(s2); m=s*INV; so=np.outer(s,s)
    rho=np.clip(G/so,-1,1); C=so*((np.sqrt(1-rho*rho)+rho*(np.pi-np.arccos(rho)))/(2*np.pi))-np.outer(m,m)
    C[np.diag_indices(n)]=s2*(0.5-1/(2*np.pi)); C=C.astype(np.float32)
    M=(G*0.5).astype(np.float32); M[np.diag_indices(n)]=0
    sources[0]=(M,(INV/s).astype(np.float32),m.astype(np.float32),
                (s2*s*(2*INV-1.5*INV+2*INV**3)).astype(np.float32),np.full(n,0.5,np.float32))
    Sig=C; m=m.astype(np.float32)
    for l in range(1,L):
        Wl=W[l]; mu=(m@Wl).astype(np.float64); Sz=((Wl.T@Sig)@Wl).astype(np.float64)
        Sz=Sz*(1.0+V[l-1])                                   # the whole calibration, one multiply
        sz2=np.diag(Sz).copy(); sz=np.sqrt(np.maximum(sz2,1e-30)); t=mu/sz
        ph=pdf(t); Ph=ndtr(t); Psi=mu*Ph+sz*ph; m2=(mu*mu+sz2)*Ph+mu*sz*ph; s3=sz2*sz; s4=sz2*sz2
        U=Wl; Ms,ps,qs,k3hs,d1s=sources[l-1]; Nt=Ms.T@U; UN=U*Nt
        k3=((ps[:,None]*np.float32(3))*UN*Nt+(qs[:,None]*np.float32(3))*UN*U+k3hs[:,None]*(U*U*U)).sum(0)
        X1=(UN*2)*ps[:,None]+(U*U)*qs[:,None]
        X2=(Nt*Nt)*ps[:,None]+(UN*2)*qs[:,None]+(U*U)*k3hs[:,None]
        K21=(X1.T@Nt+X2.T@U).astype(np.float64); np.fill_diagonal(K21,0.0)
        k3=k3.astype(np.float64)+al*A_S[l-1]*t*s3; g4=G_S[l-1]
        z=Zs[l]; a=z.mean(0); zc=z-a; k3s=(zc**3).mean(0).astype(np.float64)
        rr=k3s-k3; nz=6*np.mean(s3*s3)/N; sg=max(np.mean(rr**2)-nz,0); k3=k3+(sg/(sg+nz))*rr
        muh=mu/np.linalg.norm(mu); Ap=zc@muh.astype(np.float32)
        u=0.75*(C_S[l-1]*s3)+0.25*(((zc*zc).T@Ap).astype(np.float64)/N-K21@muh)
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
        Cn[np.diag_indices(n)]=((m2+dE2)-m_next*m_next).astype(np.float32)
        Msrc=Sz*Ph[:,None]; np.fill_diagonal(Msrc,0.0); m3=mu*m2+2*sz2*Psi
        sources[l]=(Msrc.astype(np.float32),(ph/sz).astype(np.float32),(2*Psi*(1-Ph)).astype(np.float32),
                    (m3-3*m2*Psi+2*Psi**3).astype(np.float32),Ph.astype(np.float32))
        Sig=Cn.astype(np.float32); m=m_next.astype(np.float32)
    return float(np.mean((m.astype(np.float64)-Y[-1])**2))

if __name__=="__main__":
    FIT=list(range(8)); HOLD=list(range(8,16))
    def sc(V,nets,al=1.0): return float(np.mean([predict(k,V,al=al) for k in nets]))
    print("Per-layer multiplicative calibration of Sz.  Cost: one n^2 multiply per layer, 0.0% of budget.\n")
    print(f"{'variant':34s} {'scale':>6}  {'fit(0-7)':>10} {'held(8-15)':>11}  {'vs deployed':>12}")
    base_f=sc(np.zeros(15),FIT); base_h=sc(np.zeros(15),HOLD)
    print(f"{'none (deployed)':34s} {'-':>6}  {base_f:10.4e} {base_h:11.4e}  {'1.0000':>12}")
    for name,V0 in (("cumulative  v = 2 e_l",V_CUM),("incremental v = 2 (e_l - e_{l-1})",V_INC)):
        for g in (0.5,1.0,2.0,4.0,8.0):
            f=sc(V0*g,FIT); h=sc(V0*g,HOLD)
            print(f"{name:34s} {g:6.1f}  {f:10.4e} {h:11.4e}  {h/base_h:12.4f}", flush=True)
    