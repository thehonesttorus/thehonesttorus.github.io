# Level 1 of the operator-system programme: is the propagated state positive, and does projecting it
# back into the positive cone help?  Reports min eigenvalue of the propagated post-activation covariance
# at each layer, and compares MSE with/without a PSD projection of the state.
import sys, time, numpy as np
from scipy.special import ndtr
import twin2 as T, twin4
n=1024; L=16; INV=T.INV; pdf=T.pdf
def run(W, Y, psd=0, report=False, K=3, N_MC=2400, seed=0, drop_k21=False):
    rng=np.random.default_rng(seed); X=rng.standard_normal((N_MC,n)).astype(np.float32); Zs=[]; h=X
    for l in range(L): z=h@W[l]; Zs.append(z); h=np.maximum(z,0)
    outs=[]; sources={}; info=[]
    W0=W[0]; G=W0.T@W0; s2=np.diag(G); s=np.sqrt(s2); m=s*INV; so=np.outer(s,s)
    rho=np.clip(G/so,-1,1); C=so*((np.sqrt(1-rho*rho)+rho*(np.pi-np.arccos(rho)))/(2*np.pi))-np.outer(m,m)
    C[np.diag_indices(n)]=s2*(0.5-1/(2*np.pi)); M=G*0.5; M[np.diag_indices(n)]=0
    sources[0]=(M,INV/s,m,s2*s*(2*INV-1.5*INV+2*INV**3),np.full(n,0.5)); outs.append(m); Sig=C
    for l in range(1,L):
        Wl=W[l]; mu=m@Wl; Sz=(Wl.T@Sig)@Wl; sz2=np.diag(Sz).copy(); sz=np.sqrt(sz2); t=mu/sz
        ph=pdf(t); Ph=ndtr(t); Psi=mu*Ph+sz*ph; m2=(mu*mu+sz2)*Ph+mu*sz*ph; s3=sz2*sz; s4=sz2*sz2
        k3=np.zeros(n); K21m=np.zeros((n,n)); s_lo=max(l-K,0); U=Wl.astype(np.float64)
        for src in range(l-1,s_lo-1,-1):
            Ms,ps,qs,k3hs,d1s=sources[src]; Nt=Ms.T@U; UN=U*Nt
            k3+=((ps[:,None]*3)*UN*Nt+(qs[:,None]*3)*UN*U+k3hs[:,None]*(U*U*U)).sum(0)
            K21m+=2*(UN*ps[:,None]).T@Nt+((Nt*Nt)*ps[:,None]).T@U+((U*U)*qs[:,None]).T@Nt+2*(UN*qs[:,None]).T@U+((U*U)*k3hs[:,None]).T@U
            if src>s_lo: U=W[src]@(d1s[:,None]*U)
        z=Zs[l].astype(np.float64); a=z.mean(0); zc=z-a; v=(zc*zc).mean(0)
        k3s=(zc**3).mean(0); k4s=(zc**4).mean(0)-3*v*v
        r=(k3s-k3)/s3; Xb=np.column_stack([t,np.ones(n)]); cf=np.linalg.lstsq(Xb,r,rcond=None)[0]; k3=k3+(Xb@cf)*s3
        resid=k3s-k3; noise=6*np.mean(s3*s3)/N_MC; sig=max(np.mean(resid**2)-noise,0); k3=k3+(sig/(sig+noise))*resid
        k4=np.mean(k4s/s4)*s4
        muh=mu/np.linalg.norm(mu); Ap=zc@muh; u=(zc*zc).T@Ap/N_MC-K21m@muh; K21=K21m+np.outer(u,muh)
        np.fill_diagonal(K21,0)
        t2=t*t; he2=t2-1; he4=(t2-6)*t2+3
        dE1=-k3*t*ph/(6*sz2)+(k3*k3)*he4*ph/(72*s4*sz)+k4*he2*ph/(24*s3)
        dE2=k3*ph/(3*sz)+(k3*k3)*(3*t-t*t2)*ph/(36*s4)-k4*t*ph/(12*sz2)
        m_next=Psi+dE1
        so=np.outer(sz,sz); rho=Sz/so; ds=T.hermite_d(t,ph,Ph,8); rk=rho.copy(); C=np.zeros((n,n)); fact=1.0
        for k in range(1,9):
            fact*=k; C+=rk*np.outer(ds[k-1],ds[k-1])/fact
            if k<8: rk*=rho
        C*=so
        Cg=C.copy()                                   # the Hermite (Schur-positive) part, exact diagonal below
        F2=ph/sz
        if not drop_k21: C=C+0.5*(K21*np.outer(F2,Ph)+K21.T*np.outer(Ph,F2))
        C[np.diag_indices(n)]=(m2+dE2)-m_next*m_next
        if report or psd:
            ev=np.linalg.eigvalsh(C); tr=ev.sum()
            evg=np.linalg.eigvalsh(Cg) if report else None
            if report:
                Cd=Cg.copy(); Cd[np.diag_indices(n)]=(m2+dE2)-m_next*m_next
                evd=np.linalg.eigvalsh(Cd)
                info.append((l+1,ev[0]/tr,evd[0]/tr,evg[0]/evg.sum(),np.abs(rho-np.diag(np.diag(rho))).max(),(ev<0).sum()))
        if psd:
            ev,V=np.linalg.eigh(C); floor=psd*max(ev.max(),0)
            ev=np.maximum(ev,floor); C=(V*ev)@V.T
        Msrc=Sz*Ph[:,None]; Msrc[np.diag_indices(n)]=0
        m3=mu*m2+2*sz2*Psi; sources[l]=(Msrc,ph/sz,2*Psi*(1-Ph),m3-3*m2*Psi+2*Psi**3,Ph)
        m=m_next; Sig=C; outs.append(m)
    return np.array(outs), info
if __name__=='__main__':
    seeds,Yall=T.load_nets(); W=T.regen(seeds[0]); Y=Yall[0]
    P,info=run(W,Y,report=True)
    print("Positivity of the propagated post-activation covariance (network 0, K=3).")
    print("lam_min/trace: full state | Hermite part + exact diagonal | Hermite part alone (Schur-positive) | max|rho| off-diag | #neg eigs\n")
    for l,a,b,c,mr,nn in info:
        print(f"  layer {l:2d}:  {a:+.3e}   {b:+.3e}   {c:+.3e}   {mr:.3f}   {nn}")
    print()
    for k in [0,1,2]:
        W=T.regen(seeds[k]); Y=Yall[k]
        for psd in [0, 1e-12, 1e-6, 1e-4]:
            P,_=run(W,Y,psd=psd); m=np.mean((P-Y)**2,axis=1)[-1]
            print(f"net {k} PSD projection floor={psd:g}: final MSE {m:.3e}", flush=True)
