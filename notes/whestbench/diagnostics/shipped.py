# If the reduced memory kernel's coefficient is universal in depth, ship it and drop the fit.
# Honest protocol: fit the schedule on networks 0-7, evaluate on the held-out networks 8-15.
import numpy as np, sys
from scipy.special import ndtr
import twin2 as T
n=1024; L=16; INV=T.INV; pdf=T.pdf
seeds,Yall=T.load_nets()
def run(k, mode, sched=None, N_MC=2400, K=1, seed=0, collect=False):
    """mode: 'fit' = per-network Monte-Carlo fit (current estimator);
             'ship' = shipped schedule, no Monte-Carlo at all;
             'ship+mc' = shipped schedule for the kernel, Monte-Carlo only for shrinkage and two-point."""
    W=T.regen(seeds[k]); Y=Yall[k]
    use_mc = mode!='ship'
    if use_mc:
        rng=np.random.default_rng(seed); X=rng.standard_normal((N_MC,n)).astype(np.float32); Zs=[]; h=X
        for l in range(L): z=h@W[l]; Zs.append(z); h=np.maximum(z,0)
    outs=[]; sources={}; got=[]
    W0=W[0]; G=W0.T@W0; s2=np.diag(G); s=np.sqrt(s2); m=s*INV; so=np.outer(s,s)
    rho=np.clip(G/so,-1,1); C=so*((np.sqrt(1-rho*rho)+rho*(np.pi-np.arccos(rho)))/(2*np.pi))-np.outer(m,m)
    C[np.diag_indices(n)]=s2*(0.5-1/(2*np.pi)); M=G*0.5; M[np.diag_indices(n)]=0
    sources[0]=(M,INV/s,m,s2*s*(2*INV-1.5*INV+2*INV**3),np.full(n,0.5)); outs.append(m); Sig=C
    for l in range(1,L):
        Wl=W[l]; mu=m@Wl; Sz=(Wl.T@Sig)@Wl; sz2=np.diag(Sz).copy(); sz=np.sqrt(sz2); t=mu/sz
        ph=pdf(t); Ph=ndtr(t); Psi=mu*Ph+sz*ph; m2=(mu*mu+sz2)*Ph+mu*sz*ph; s3=sz2*sz; s4=sz2*sz2
        k3=np.zeros(n); U=Wl.astype(np.float64); s_lo=max(l-K,0)
        for src in range(l-1,s_lo-1,-1):
            Ms,ps,qs,k3hs,d1s=sources[src]; Nt=Ms.T@U; UN=U*Nt
            k3+=((ps[:,None]*3)*UN*Nt+(qs[:,None]*3)*UN*U+k3hs[:,None]*(U*U*U)).sum(0)
            if src>s_lo: U=W[src]@(d1s[:,None]*U)
        if use_mc:
            z=Zs[l].astype(np.float64); a=z.mean(0); zc=z-a; v=(zc*zc).mean(0)
            k3s=(zc**3).mean(0); k4s=(zc**4).mean(0)-3*v*v
        if mode=='fit':
            r=(k3s-k3)/s3; Xb=np.column_stack([t,np.ones(n)]); cf=np.linalg.lstsq(Xb,r,rcond=None)[0]
            g4=np.mean(k4s/s4); k3=k3+(Xb@cf)*s3
            if collect: got.append((cf[0],g4))
        else:
            av,g4=sched[l-1]; k3=k3+av*t*s3
        if use_mc:
            resid=k3s-k3; noise=6*np.mean(s3*s3)/N_MC; sig=max(np.mean(resid**2)-noise,0)
            k3=k3+(sig/(sig+noise))*resid
        k4=g4*s4
        t2=t*t; he2=t2-1; he4=(t2-6)*t2+3
        dE1=-k3*t*ph/(6*sz2)+(k3*k3)*he4*ph/(72*s4*sz)+k4*he2*ph/(24*s3)
        dE2=k3*ph/(3*sz)+(k3*k3)*(3*t-t*t2)*ph/(36*s4)-k4*t*ph/(12*sz2)
        m_next=Psi+dE1
        so=np.outer(sz,sz); rho=Sz/so; ds=T.hermite_d(t,ph,Ph,8); rk=rho.copy(); C=np.zeros((n,n)); fact=1.0
        for kk in range(1,9):
            fact*=kk; C+=rk*np.outer(ds[kk-1],ds[kk-1])/fact
            if kk<8: rk*=rho
        C*=so
        if use_mc:
            muh=mu/np.linalg.norm(mu); Ap=zc@muh; u=(zc*zc).T@Ap/N_MC
            K21=np.outer(u,muh); np.fill_diagonal(K21,0); F2=ph/sz
            C=C+0.5*(K21*np.outer(F2,Ph)+K21.T*np.outer(Ph,F2))
        C[np.diag_indices(n)]=(m2+dE2)-m_next*m_next
        Msrc=Sz*Ph[:,None]; Msrc[np.diag_indices(n)]=0
        m3=mu*m2+2*sz2*Psi; sources[l]=(Msrc,ph/sz,2*Psi*(1-Ph),m3-3*m2*Psi+2*Psi**3,Ph)
        m=m_next; Sig=C; outs.append(m)
    P=np.array(outs)
    return (np.mean((P-Y)**2,axis=1)[-1], got)
# --- fit the schedule on networks 0-7, averaging over three Monte-Carlo seeds to beat down noise ---
acc=[]
for k in range(8):
    for sd in [0,1,2]:
        _,g=run(k,'fit',seed=sd,collect=True); acc.append(g)
acc=np.array(acc)                     # (samples, layer, 2)
sched=acc.mean(0)
print("Shipped schedule fitted on networks 0-7 (24 fits): coefficient of t*sigma^3, and mean excess kurtosis")
for l in range(L-1): print(f"  layer {l+2:2d}:  a = {sched[l,0]:+.5f}   g4 = {sched[l,1]:+.5f}")
print("\nHeld-out networks 8-15:")
print(" net |  MC-fitted kernel  |  shipped, no MC   |  shipped + MC channels")
r1=[];r2=[];r3=[]
for k in range(8,16):
    a,_=run(k,'fit'); b,_=run(k,'ship',sched=sched); c,_=run(k,'ship+mc',sched=sched)
    r1.append(a); r2.append(b); r3.append(c)
    print(f" {k:3d} |     {a:.3e}      |     {b:.3e}     |     {c:.3e}", flush=True)
print(f"mean |     {np.mean(r1):.3e}      |     {np.mean(r2):.3e}     |     {np.mean(r3):.3e}")
