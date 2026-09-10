# OQCP memory channel: carry the long-lived cubic sector as R symmetric rank-one modes
#   T ~ sum_r lambda_r v_r^{(x)3},   readout  k3_mem,j = sum_r lambda_r (W^T v_r)_j^3
# Transport is one matvec per mode per layer and is shared with the readout:
#   U = W_l^T V   (used for readout),   V <- Phi_l * U   (transport to the next layer)
# Sources of age 1 stay exact (the existing channel); the memory carries age >= 2 only.
import numpy as np
from scipy.special import ndtr
import twin2 as T
n=1024; L=16; INV=T.INV; pdf=T.pdf
seeds,Yall=T.load_nets()
SCHED=np.array([[0.00094,0.00715],[0.00783,0.01340],[0.01134,0.01947],[0.01406,0.02302],[0.01658,0.02679],
 [0.01853,0.03065],[0.02040,0.03478],[0.02217,0.03824],[0.02331,0.04162],[0.02464,0.04615],
 [0.02555,0.04940],[0.02684,0.05096],[0.02781,0.05280],[0.02838,0.05462],[0.02940,0.05873]])
def run(k, R=0, B=256, kernel=True, mc=True, K=1, N_MC=2400, seed=0, polar=False, sel='mag'):
    W=T.regen(seeds[k]); Y=Yall[k]
    if mc:
        rng=np.random.default_rng(seed); X=rng.standard_normal((N_MC,n)).astype(np.float32); Zs=[]; h=X
        for l in range(L): z=h@W[l]; Zs.append(z); h=np.maximum(z,0)
    outs=[]; sources={}; V=np.zeros((n,0)); lam=np.zeros(0); pending=None
    W0=W[0]; G=W0.T@W0; s2=np.diag(G); s=np.sqrt(s2); m=s*INV; so=np.outer(s,s)
    rho=np.clip(G/so,-1,1); C=so*((np.sqrt(1-rho*rho)+rho*(np.pi-np.arccos(rho)))/(2*np.pi))-np.outer(m,m)
    C[np.diag_indices(n)]=s2*(0.5-1/(2*np.pi)); M=G*0.5; M[np.diag_indices(n)]=0
    k3h0=s2*s*(2*INV-1.5*INV+2*INV**3)
    sources[0]=(M,INV/s,m,k3h0,np.full(n,0.5)); outs.append(m); Sig=C
    for l in range(1,L):
        Wl=W[l]; mu=m@Wl; Sz=(Wl.T@Sig)@Wl; sz2=np.diag(Sz).copy(); sz=np.sqrt(sz2); t=mu/sz
        ph=pdf(t); Ph=ndtr(t); Psi=mu*Ph+sz*ph; m2=(mu*mu+sz2)*Ph+mu*sz*ph; s3=sz2*sz; s4=sz2*sz2
        k3=np.zeros(n)
        if K>0:                                     # exact age-1 source
            U0=Wl.astype(np.float64); s_lo=max(l-K,0)
            for src in range(l-1,s_lo-1,-1):
                Ms,ps,qs,k3hs,d1s=sources[src]; Nt=Ms.T@U0; UN=U0*Nt
                k3+=((ps[:,None]*3)*UN*Nt+(qs[:,None]*3)*UN*U0+k3hs[:,None]*(U0*U0*U0)).sum(0)
                if src>s_lo: U0=W[src]@(d1s[:,None]*U0)
        # ---- memory channel: readout, then transport, then absorb the source that just aged out ----
        Um=None
        if R>0:
            if V.shape[1]>0:
                Um=Wl.T@V                            # (n outputs, R modes)
                k3+=(Um**3)@lam
            if kernel: k3=k3+SCHED[l-1,0]*t*s3*0.0   # kernel replaced by the explicit memory
        elif kernel:
            k3=k3+SCHED[l-1,0]*t*s3
        g4=SCHED[l-1,1]
        if mc:
            z=Zs[l].astype(np.float64); a=z.mean(0); zc=z-a; v=(zc*zc).mean(0)
            k3s=(zc**3).mean(0); k4s=(zc**4).mean(0)-3*v*v
            resid=k3s-k3; noise=6*np.mean(s3*s3)/N_MC; sig=max(np.mean(resid**2)-noise,0)
            k3=k3+(sig/(sig+noise))*resid
        k4v=g4*s4
        t2=t*t; he2=t2-1; he4=(t2-6)*t2+3
        dE1=-k3*t*ph/(6*sz2)+(k3*k3)*he4*ph/(72*s4*sz)+k4v*he2*ph/(24*s3)
        dE2=k3*ph/(3*sz)+(k3*k3)*(3*t-t*t2)*ph/(36*s4)-k4v*t*ph/(12*sz2)
        m_next=Psi+dE1
        so=np.outer(sz,sz); rhoc=Sz/so; ds=T.hermite_d(t,ph,Ph,8); rk=rhoc.copy(); C=np.zeros((n,n)); fact=1.0
        for kk in range(1,9):
            fact*=kk; C+=rk*np.outer(ds[kk-1],ds[kk-1])/fact
            if kk<8: rk*=rhoc
        C*=so
        if mc:
            muh=mu/np.linalg.norm(mu); Ap=zc@muh; u=(zc*zc).T@Ap/N_MC
            K21=np.outer(u,muh); np.fill_diagonal(K21,0); F2=ph/sz
            C=C+0.5*(K21*np.outer(F2,Ph)+K21.T*np.outer(Ph,F2))
        C[np.diag_indices(n)]=(m2+dE2)-m_next*m_next
        Msrc=Sz*Ph[:,None]; Msrc[np.diag_indices(n)]=0
        m3=mu*m2+2*sz2*Psi; k3h=m3-3*m2*Psi+2*Psi**3
        sources[l]=(Msrc,ph/sz,2*Psi*(1-Ph),k3h,Ph)
        if R>0:
            if V.shape[1]>0: V=Ph[:,None]*Um          # transport existing modes into h_{l+1} space
            if pending is not None:                   # absorb the source that has now aged past 1
                pv,pl=pending
                V=np.concatenate([V,pv],1) if V.shape[1] else pv
                lam=np.concatenate([lam,pl])
            # this layer's own source becomes pending.  The source has three diagrams; only the
            # diagonal one is natively lambda*v^(x)3, so polarize the other two:
            #   3 sym(e,b,b) = [(e+b)^3 + (e-b)^3]/2 - e^3 ,  3 sym(e,e,b) = [(e+b)^3 - (e-b)^3]/2 - b^3
            # with b = m_beta / ||m_beta||, m_beta the beta-th column of the source matrix Msrc.
            pv_=ph/sz; qv_=2*Psi*(1-Ph)
            if polar:
                Mb=Msrc.copy(); nb=np.linalg.norm(Mb,axis=0); nb[nb==0]=1.0; Bm=Mb/nb
                lam_e = k3h - pv_*nb*nb
                lam_b = -qv_*nb
                lam_p = (pv_*nb*nb + qv_*nb)/2.0
                lam_m = (pv_*nb*nb - qv_*nb)/2.0
                E=np.eye(n)
                cand_v=np.concatenate([E, Bm, (E+Bm)/np.sqrt(2.0), (E-Bm)/np.sqrt(2.0)],1)
                r2=2.0*np.sqrt(2.0)
                cand_l=np.concatenate([lam_e, lam_b, lam_p*r2, lam_m*r2])
            else:
                cand_v=np.eye(n); cand_l=k3h.copy()
            idx=np.argsort(-np.abs(cand_l))[:B]
            pending=(cand_v[:,idx].copy(), cand_l[idx].copy())
            if V.shape[1]>R:                          # compress back to R modes
                if sel=='mag':
                    w=np.abs(lam)*np.linalg.norm(V,axis=0)**3
                    keep=np.argsort(-w)[:R]; V=V[:,keep]; lam=lam[keep]
        m=m_next; Sig=C; outs.append(m)
    return np.mean((np.array(outs)-Y)**2,axis=1)[-1]
if __name__=="__main__":
    import sys
    print("Held-out networks 8-15, final-layer MSE.\n")
    print(f"{'config':46s}   mean         per network")
    for name,kw in [("current: age-1 exact + shipped kernel",  dict(R=0,kernel=True)),
                    ("polarized memory R=128, B=1024",         dict(R=128,B=1024,polar=True)),
                    ("polarized memory R=256, B=2048",         dict(R=256,B=2048,polar=True)),
                    ("polarized memory R=512, B=4096",         dict(R=512,B=4096,polar=True)),
                    ("polarized memory R=1024, B=4096",        dict(R=1024,B=4096,polar=True))]:
        v=[run(k,**kw) for k in range(8,16)]
        print(f"{name:46s} {np.mean(v):.3e}   " + " ".join(f"{x:.2e}" for x in v), flush=True)
