# Per-layer single-step readout error with EXACT means and PROPER cumulant ratios from N=1e6
# (lambda_3 = kappa_3/sigma^3 built from raw moments, not from a sample-standardised Hermite average).
# Layer 1 is the control: there the covariance is the closed-form arcsine matrix, so whatever the
# error is at layer 1 is the readout, and everything above it is the propagated covariance.
import numpy as np, sys, os
from scipy.special import ndtr
import ablate as AB, twin2 as T
n=1024; L=16; INV=T.INV; pdf=T.pdf
seeds,Yall=AB.seeds,AB.Yall
A_S,G_S,C_S=AB.A_S,AB.G_S,AB.C_S
NB=1000000; NETS=[0,1]
def run(k,lam,N=2400,seed=0):
    W=AB.getW(k); Y=Yall[k]
    rng=np.random.default_rng(seed); X=rng.standard_normal((N,n)).astype(np.float32); Zs=[]; h=X
    for l in range(L): z=h@W[l]; Zs.append(z); h=np.maximum(z,0)
    sources={}; prof=[]
    W0=W[0].astype(np.float64); G=W0.T@W0; s2=np.diag(G); s=np.sqrt(s2); m=s*INV; so=np.outer(s,s)
    rho=np.clip(G/so,-1,1); C=so*((np.sqrt(1-rho*rho)+rho*(np.pi-np.arccos(rho)))/(2*np.pi))-np.outer(m,m)
    C[np.diag_indices(n)]=s2*(0.5-1/(2*np.pi))
    M=(G*0.5).astype(np.float32); M[np.diag_indices(n)]=0
    sources[0]=(M,(INV/s).astype(np.float32),m.astype(np.float32),
                (s2*s*(2*INV-1.5*INV+2*INV**3)).astype(np.float32),np.full(n,0.5,np.float32))
    Sig=C.astype(np.float32); m=Y[0].astype(np.float32)
    for l in range(1,L):
        Wl=W[l]; mu=(Y[l-1]@Wl.astype(np.float64))
        Sz=((Wl.T@Sig)@Wl).astype(np.float64)
        sz2=np.diag(Sz).copy(); sz=np.sqrt(np.maximum(sz2,1e-30)); t=mu/sz
        ph=pdf(t); Ph=ndtr(t); Psi=mu*Ph+sz*ph; m2=(mu*mu+sz2)*Ph+mu*sz*ph; s3=sz2*sz; s4=sz2*sz2
        U=Wl; Ms,ps,qs,k3hs,d1s=sources[l-1]; Nt=Ms.T@U; UN=U*Nt
        k3m=((ps[:,None]*np.float32(3))*UN*Nt+(qs[:,None]*np.float32(3))*UN*U+k3hs[:,None]*(U*U*U)).sum(0)
        X1=(UN*2)*ps[:,None]+(U*U)*qs[:,None]
        X2=(Nt*Nt)*ps[:,None]+(UN*2)*qs[:,None]+(U*U)*k3hs[:,None]
        K21=(X1.T@Nt+X2.T@U).astype(np.float64); np.fill_diagonal(K21,0.0)
        k3m=k3m.astype(np.float64)+A_S[l-1]*t*s3
        z=Zs[l]; a=z.mean(0); zc=z-a; k3s=(zc**3).mean(0).astype(np.float64)
        rr=k3s-k3m; nz=6*np.mean(s3*s3)/N; sg=max(np.mean(rr**2)-nz,0); k3m=k3m+(sg/(sg+nz))*rr
        k3=lam[0][l]*s3; k4v=lam[1][l]*s4              # oracle cumulants
        t2=t*t; he2=t2-1; he4=(t2-6)*t2+3
        def dE(K3,K4): return -K3*t*ph/(6*sz2)+(K3*K3)*he4*ph/(72*s4*sz)+K4*he2*ph/(24*s3)
        R_or=Y[l]-(Psi+dE(k3,k4v))                      # exact cumulants
        R_md=Y[l]-(Psi+dE(k3m,G_S[l-1]*s4))             # our own cumulants
        R_ga=Y[l]-Psi                                   # no cumulant correction at all
        prof.append((np.mean(R_ga**2),np.mean(R_md**2),np.mean(R_or**2),
                     float(np.mean((lam[2][l]/sz-1)))))
        # continue with the DEPLOYED model so Sigma is the propagated one
        muh=mu/np.linalg.norm(mu); Ap=zc@muh.astype(np.float32)
        u=0.75*(C_S[l-1]*s3)+0.25*(((zc*zc).T@Ap).astype(np.float64)/N-K21@muh)
        K21=K21+np.outer(u,muh); np.fill_diagonal(K21,0.0)
        k4m=G_S[l-1]*s4
        dE1=dE(k3m,k4m)
        dE2=k3m*ph/(3*sz)+(k3m*k3m)*(3*t-t*t2)*ph/(36*s4)-k4m*t*ph/(12*sz2)
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
        Sig=Cn.astype(np.float32)
    return prof
LAM={}
for k in NETS:
    d=np.load(f"mom_{k}_{NB}.npz"); LAM[k]=(d['lam3'],d['lam4'],d['sd'])
P=[run(k,LAM[k]) for k in NETS]
print("Single-step readout MSE with EXACT layer means (networks 0,1; oracle cumulants from 1e6 samples).")
print("Layer 1 has the exact arcsine covariance; everything above it carries the propagated one.\n")
print(f"{'layer':>5} {'Gaussian':>11} {'our cumulants':>14} {'exact cumulants':>16} "
      f"{'cum. share':>11} {'sigma err':>10}")
for i in range(L-1):
    g=np.mean([p[i][0] for p in P]); md=np.mean([p[i][1] for p in P]); orc=np.mean([p[i][2] for p in P])
    se=np.mean([p[i][3] for p in P])
    print(f"{i+1:5d} {g:11.3e} {md:14.3e} {orc:16.3e} {1-orc/md:11.3f} {se:10.2e}")
