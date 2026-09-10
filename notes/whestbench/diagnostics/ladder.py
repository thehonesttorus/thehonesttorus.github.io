# plain K2 vs renormalized K2 vs factorized K3 vs low-rank memory vs full.
# Schedule fitted on networks 0-7 (from shipped.py); evaluated on held-out networks 8-15.
import numpy as np
from scipy.special import ndtr
import twin2 as T
n=1024; L=16; INV=T.INV; pdf=T.pdf
seeds,Yall=T.load_nets()
SCHED=np.array([[0.00094,0.00715],[0.00783,0.01340],[0.01134,0.01947],[0.01406,0.02302],[0.01658,0.02679],
 [0.01853,0.03065],[0.02040,0.03478],[0.02217,0.03824],[0.02331,0.04162],[0.02464,0.04615],
 [0.02555,0.04940],[0.02684,0.05096],[0.02781,0.05280],[0.02838,0.05462],[0.02940,0.05873]])
def run(k, K=0, kernel=False, mc=False, k4=False, N_MC=2400, seed=0):
    W=T.regen(seeds[k]); Y=Yall[k]
    if mc:
        rng=np.random.default_rng(seed); X=rng.standard_normal((N_MC,n)).astype(np.float32); Zs=[]; h=X
        for l in range(L): z=h@W[l]; Zs.append(z); h=np.maximum(z,0)
    outs=[]; sources={}
    W0=W[0]; G=W0.T@W0; s2=np.diag(G); s=np.sqrt(s2); m=s*INV; so=np.outer(s,s)
    rho=np.clip(G/so,-1,1); C=so*((np.sqrt(1-rho*rho)+rho*(np.pi-np.arccos(rho)))/(2*np.pi))-np.outer(m,m)
    C[np.diag_indices(n)]=s2*(0.5-1/(2*np.pi)); M=G*0.5; M[np.diag_indices(n)]=0
    sources[0]=(M,INV/s,m,s2*s*(2*INV-1.5*INV+2*INV**3),np.full(n,0.5)); outs.append(m); Sig=C
    for l in range(1,L):
        Wl=W[l]; mu=m@Wl; Sz=(Wl.T@Sig)@Wl; sz2=np.diag(Sz).copy(); sz=np.sqrt(sz2); t=mu/sz
        ph=pdf(t); Ph=ndtr(t); Psi=mu*Ph+sz*ph; m2=(mu*mu+sz2)*Ph+mu*sz*ph; s3=sz2*sz; s4=sz2*sz2
        k3=np.zeros(n)
        if K>0:
            U=Wl.astype(np.float64); s_lo=max(l-K,0)
            for src in range(l-1,s_lo-1,-1):
                Ms,ps,qs,k3hs,d1s=sources[src]; Nt=Ms.T@U; UN=U*Nt
                k3+=((ps[:,None]*3)*UN*Nt+(qs[:,None]*3)*UN*U+k3hs[:,None]*(U*U*U)).sum(0)
                if src>s_lo: U=W[src]@(d1s[:,None]*U)
        if kernel: k3=k3+SCHED[l-1,0]*t*s3
        g4=SCHED[l-1,1] if (kernel or k4) else 0.0
        if mc:
            z=Zs[l].astype(np.float64); a=z.mean(0); zc=z-a; v=(zc*zc).mean(0)
            k3s=(zc**3).mean(0); k4s=(zc**4).mean(0)-3*v*v
            if not kernel:
                r=(k3s-k3)/s3; Xb=np.column_stack([t,np.ones(n)]); cf=np.linalg.lstsq(Xb,r,rcond=None)[0]
                k3=k3+(Xb@cf)*s3; g4=np.mean(k4s/s4)
            resid=k3s-k3; noise=6*np.mean(s3*s3)/N_MC; sig=max(np.mean(resid**2)-noise,0)
            k3=k3+(sig/(sig+noise))*resid
        k4v=g4*s4
        t2=t*t; he2=t2-1; he4=(t2-6)*t2+3
        dE1=-k3*t*ph/(6*sz2)+(k3*k3)*he4*ph/(72*s4*sz)+k4v*he2*ph/(24*s3)
        dE2=k3*ph/(3*sz)+(k3*k3)*(3*t-t*t2)*ph/(36*s4)-k4v*t*ph/(12*sz2)
        m_next=Psi+dE1
        so=np.outer(sz,sz); rho=Sz/so; ds=T.hermite_d(t,ph,Ph,8); rk=rho.copy(); C=np.zeros((n,n)); fact=1.0
        for kk in range(1,9):
            fact*=kk; C+=rk*np.outer(ds[kk-1],ds[kk-1])/fact
            if kk<8: rk*=rho
        C*=so
        if mc:
            muh=mu/np.linalg.norm(mu); Ap=zc@muh; u=(zc*zc).T@Ap/N_MC
            K21=np.outer(u,muh); np.fill_diagonal(K21,0); F2=ph/sz
            C=C+0.5*(K21*np.outer(F2,Ph)+K21.T*np.outer(Ph,F2))
        C[np.diag_indices(n)]=(m2+dE2)-m_next*m_next
        Msrc=Sz*Ph[:,None]; Msrc[np.diag_indices(n)]=0
        m3=mu*m2+2*sz2*Psi; sources[l]=(Msrc,ph/sz,2*Psi*(1-Ph),m3-3*m2*Psi+2*Psi**3,Ph)
        m=m_next; Sig=C; outs.append(m)
    return np.mean((np.array(outs)-Y)**2,axis=1)[-1]
ladder=[("plain K2 (Gaussian closure)",           dict(K=0,kernel=False,mc=False)),
        ("K2 + rank-1 memory kernel (shipped)",   dict(K=0,kernel=True, mc=False)),
        ("factorized K3, depth 1",                dict(K=1,kernel=False,mc=False)),
        ("factorized K3, depth 3",                dict(K=3,kernel=False,mc=False)),
        ("K3 depth 1 + rank-1 memory kernel",     dict(K=1,kernel=True, mc=False)),
        ("K3 depth 3 + rank-1 memory kernel",     dict(K=3,kernel=True, mc=False)),
        ("full: K3 depth 1 + MC residual + 2pt",  dict(K=1,kernel=False,mc=True))]
print("Held-out networks 8-15, final-layer MSE.  Schedule fitted on networks 0-7.\n")
print(f"{'method':40s}   mean        per network")
for name,kw in ladder:
    v=[run(k,**kw) for k in range(8,16)]
    print(f"{name:40s} {np.mean(v):.3e}   " + " ".join(f"{x:.2e}" for x in v), flush=True)
