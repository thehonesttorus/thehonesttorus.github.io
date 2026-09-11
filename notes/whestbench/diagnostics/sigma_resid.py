# With exact means AND exact kappa_3, kappa_4 the closure still leaves 9.87e-08, and the Gram-Charlier
# truncation accounts for only 1.2e-09 of that.  What is left can only be the propagated variance.
# Measure it WITHOUT Monte-Carlo noise: dPsi/dsigma = phi(t) exactly, so given the true layer mean Y_l
# and the true mu_l = Y_{l-1} W_l we can SOLVE for the sigma that the closure would need,
#     Psi(mu,sigma) + dE1(mu,sigma,k3,k4) = Y_l,
# and compare it to the sigma the closure propagates.  Then ask whether the discrepancy is a universal
# function of relative depth -- if it is, it is one shipped number per layer and it is free.
import numpy as np, sys, os
from scipy.special import ndtr
import ablate as AB, twin2 as T
n=1024; L=16; INV=T.INV; pdf=T.pdf
seeds,Yall=AB.seeds,AB.Yall
A_S,G_S,C_S=AB.A_S,AB.G_S,AB.C_S
NB=120000
NETS=[int(x) for x in sys.argv[1].split(',')] if len(sys.argv)>1 else [0,1]
USE_EXACT_CUM = os.path.exists(f"rw_{NETS[0]}_{NB}.npz")

def trace(k,N=2400,seed=0,lam=None,exact_mean=True):
    W=AB.getW(k); Y=Yall[k]
    rng=np.random.default_rng(seed); X=rng.standard_normal((N,n)).astype(np.float32); Zs=[]; h=X
    for l in range(L): z=h@W[l]; Zs.append(z); h=np.maximum(z,0)
    sources={}; rec={}
    W0=W[0]; G=W0.T@W0; s2=np.diag(G); s=np.sqrt(s2); m=s*INV; so=np.outer(s,s)
    rho=np.clip(G/so,-1,1); C=so*((np.sqrt(1-rho*rho)+rho*(np.pi-np.arccos(rho)))/(2*np.pi))-np.outer(m,m)
    C[np.diag_indices(n)]=s2*(0.5-1/(2*np.pi)); C=C.astype(np.float32)
    M=(G*0.5).astype(np.float32); M[np.diag_indices(n)]=0
    sources[0]=(M,(INV/s).astype(np.float32),m.astype(np.float32),
                (s2*s*(2*INV-1.5*INV+2*INV**3)).astype(np.float32),np.full(n,0.5,np.float32))
    Sig=C; m=(Y[0] if exact_mean else m).astype(np.float32)
    for l in range(1,L):
        Wl=W[l]; mu=(m@Wl).astype(np.float64); Sz=((Wl.T@Sig)@Wl).astype(np.float64)
        sz2=np.diag(Sz).copy(); sz=np.sqrt(np.maximum(sz2,1e-30)); t=mu/sz
        ph=pdf(t); Ph=ndtr(t); Psi=mu*Ph+sz*ph; m2=(mu*mu+sz2)*Ph+mu*sz*ph; s3=sz2*sz; s4=sz2*sz2
        U=Wl; Ms,ps,qs,k3hs,d1s=sources[l-1]; Nt=Ms.T@U; UN=U*Nt
        k3=((ps[:,None]*np.float32(3))*UN*Nt+(qs[:,None]*np.float32(3))*UN*U+k3hs[:,None]*(U*U*U)).sum(0)
        X1=(UN*2)*ps[:,None]+(U*U)*qs[:,None]
        X2=(Nt*Nt)*ps[:,None]+(UN*2)*qs[:,None]+(U*U)*k3hs[:,None]
        K21=X1.T@Nt+X2.T@U; np.fill_diagonal(K21,0.0)
        k3=k3.astype(np.float64)+A_S[l-1]*t*s3; g4=G_S[l-1]
        z=Zs[l]; a=z.mean(0); zc=z-a; k3s=(zc**3).mean(0).astype(np.float64)
        resid=k3s-k3; noise=6*np.mean(s3*s3)/N; sg=max(np.mean(resid**2)-noise,0)
        k3=k3+(sg/(sg+noise))*resid
        k4v=g4*s4
        if lam is not None: k3=lam[3][l]*s3; k4v=lam[4][l]*s4
        muh=mu/np.linalg.norm(mu); Ap=zc@muh.astype(np.float32)
        u=0.75*(C_S[l-1]*s3)+0.25*(((zc*zc).T@Ap).astype(np.float64)/N-K21@muh)
        K21=K21.astype(np.float64)+np.outer(u,muh); np.fill_diagonal(K21,0.0)
        t2=t*t; he2=t2-1; he4=(t2-6)*t2+3
        def readout(sig):
            tt=mu/sig; pp=pdf(tt); PP=ndtr(tt); S3=sig**3; S4=sig**4
            PS=mu*PP+sig*pp; T2=tt*tt
            d1=-k3*tt*pp/(6*sig*sig)+(k3*k3)*((T2-6)*T2+3)*pp/(72*S4*sig)+k4v*(T2-1)*pp/(24*S3)
            return PS+d1, pp
        # Newton on sigma:  dPsi/dsigma = phi(t)
        _val,_pp=readout(sz)
        rec[l]=dict(t=t.copy(),sz=sz.copy(),resid=(Y[l]-_val).copy(),Ph=Ph.copy())
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
        Sig=Cn.astype(np.float32)
        m=(Y[l] if exact_mean else m_next).astype(np.float32)
    return rec

LAM={}
for k in NETS:
    fn=f"rw_{k}_{NB}.npz"
    LAM[k]=({3:np.load(fn)['Hm'][:,3],4:np.load(fn)['Hm'][:,4]} if os.path.exists(fn) else None)
R={k:trace(k,lam=LAM[k]) for k in NETS}
print("Single-step residual with EXACT means and EXACT kappa_3, kappa_4.")
print("  R_l = Y_l - Psi(mu_exact, sigma_model) - dE1.   rms(R_l)^2 is that layer's readout MSE.")
print("  'dsig/sig' = R/(phi(t) sigma), the sigma error R would imply, on |t|<2.\n")
print(f"{'layer':>5} {'rms(R_l)':>10} {'MSE_l':>10} {'R2 sig-scale':>13} {'R2 +t':>8} {'R2 +t,t2':>9}"
      f" {'R2 free g(t)':>13} {'dsig/sig med':>13} {'dsig/sig sd':>12}")
for l in range(1,L):
    rows=[]
    for k in NETS:
        d=R[k][l]; t=d['t']; sz=d['sz']; ph=pdf(t)
        Y=Yall[k][l]; Ph=ndtr(t); mu=t*sz
        Psi=mu*Ph+sz*ph
        Rl=d['resid']
        b1=(ph*sz)[:,None]
        b2=np.column_stack([ph*sz,ph*sz*t]); b3=np.column_stack([ph*sz,ph*sz*t,ph*sz*t*t])
        # free g(t): 12 bins in t, i.e. the best possible correction that is a function of t and sigma
        qs=np.quantile(t,np.linspace(0,1,13)); idx=np.clip(np.searchsorted(qs,t)-1,0,11)
        B=np.zeros((n,12)); B[np.arange(n),idx]=ph*sz
        def r2(X):
            c,_,_,_=np.linalg.lstsq(X,Rl,rcond=None); r=Rl-X@c; return 1-np.sum(r*r)/np.sum(Rl*Rl)
        msk=np.abs(t)<2.0
        ds=Rl[msk]/(ph[msk]*sz[msk])/sz[msk]
        rows.append([np.sqrt(np.mean(Rl*Rl)),np.mean(Rl*Rl),r2(b1),r2(b2),r2(b3),r2(B),
                     np.median(ds),np.std(ds)])
    m=np.mean(rows,0)
    print(f"{l:5d} {m[0]:10.3e} {m[1]:10.3e} {m[2]:13.4f} {m[3]:8.4f} {m[4]:9.4f} {m[5]:13.4f}"
          f" {m[6]:13.2e} {m[7]:12.2e}")
