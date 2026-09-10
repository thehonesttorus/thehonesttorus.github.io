# v4 twin: + low-rank two-point third-cumulant correction of the off-diagonal covariance closure
#   K21 (kappa(z_a,z_a,z_b)) = recent-source model (K sources) + rank-r SVD of the MC residual
import sys, time, numpy as np
from scipy.special import ndtr
import twin2 as T
n=1024; L=16; INV=T.INV; pdf=T.pdf
def predict(W, N_MC=4000, K=3, T_ACT=3.0, KH=8, seed=0, rank=16, k21_model=True, k21_mc=True, k21_mu=False, k3_mu=False, trank=0, K0=2, K21_K=99, site_frac=1.0, site_pow=3.0, shrink_k3='auto', k4mode='mean', var_corr=True, oracle_cov=None, oracle_mean=None, k3_fit=True):
    rng=np.random.default_rng(seed); Zs=None
    if N_MC>0:
        X=rng.standard_normal((N_MC,n)).astype(np.float32); Zs=[]; h=X
        for l in range(L): z=h@W[l]; Zs.append(z); h=np.maximum(z,0)
    outs=[]; sources={}
    W0=W[0]; G=W0.T@W0; s2=np.diag(G); s=np.sqrt(s2); m=s*INV; so=np.outer(s,s)
    rho=np.clip(G/so,-1,1); C=so*((np.sqrt(1-rho*rho)+rho*(np.pi-np.arccos(rho)))/(2*np.pi))-np.outer(m,m)
    C[np.diag_indices(n)]=s2*(0.5-1/(2*np.pi)); M=G*0.5; M[np.diag_indices(n)]=0
    sources[0]=(M,INV/s,m,s2*s*(2*INV-1.5*INV+2*INV**3),np.full(n,0.5,dtype=np.float32)); outs.append(m); Sig=C
    for l in range(1,L):
        if oracle_cov is not None: Sig=oracle_cov[l]
        if oracle_mean is not None: m=oracle_mean[l-1]
        Wl=W[l]; mu=m@Wl; Sz=(Wl.T@Sig)@Wl; sz2=np.diag(Sz).copy(); sz=np.sqrt(sz2); t=mu/sz; ph=pdf(t); Ph=ndtr(t)
        Psi=mu*Ph+sz*ph; m2=(mu*mu+sz2)*Ph+mu*sz*ph; s3=sz2*sz; s4=sz2*sz2
        # ---- recent-source model: per-neuron k3 and two-point K21 ----
        k3=np.zeros(n); K21m=np.zeros((n,n))
        if K>0:
            s_lo=max(l-K,0); U=Wl.astype(np.float64); age=0
            for src in range(l-1,s_lo-1,-1):
                if trank>0 and age>=K0:
                    Uu,sv,Vt=np.linalg.svd(U,full_matrices=False); U=(Uu[:,:trank]*sv[:trank])@Vt[:trank]
                Ms,ps,qs,k3hs,d1s=sources[src]
                if site_frac<1.0 and age>=K0:
                    # Cartan-compatible truncation: restrict to a sub-unit-space Y (corner 1_Y A 1_Y).
                    imp=np.abs(ps)*np.linalg.norm(U,axis=1)**site_pow
                    q=max(1,int(site_frac*n)); keep=np.argpartition(-imp,q-1)[:q]
                    mask=np.zeros(n); mask[keep]=1.0; U=U*mask[:,None]
                Nt=Ms.T@U; UN=U*Nt
                k3+=((ps[:,None]*3)*UN*Nt+(qs[:,None]*3)*UN*U+k3hs[:,None]*(U*U*U)).sum(0)
                if k21_model and age<K21_K:
                    K21m+=2*(UN*ps[:,None]).T@Nt+((Nt*Nt)*ps[:,None]).T@U+((U*U)*qs[:,None]).T@Nt+2*(UN*qs[:,None]).T@U+((U*U)*k3hs[:,None]).T@U
                if src>s_lo: U=W[src]@(d1s[:,None]*U); age+=1
        k4=np.zeros(n); K21=K21m
        if Zs is not None:
            zl=Zs[l].astype(np.float64); a=zl.mean(0); zc=zl-a; v=(zc*zc).mean(0); k3s=(zc**3).mean(0); k4s=(zc**4).mean(0)-3*v*v
            if k3_fit:
                r=(k3s-k3)/s3; Xb=np.column_stack([t,np.ones(n)]); cf=np.linalg.lstsq(Xb,r,rcond=None)[0]; k3=k3+(Xb@cf)*s3
                if shrink_k3=='auto':
                    resid=k3s-k3; noise=6*s3*s3/N_MC; sig=max(np.mean(resid**2)-np.mean(noise),0); w=sig/(sig+np.mean(noise)); k3=k3+w*resid
                elif shrink_k3>0: k3=k3+shrink_k3*(k3s-k3)
            g4=k4s/s4
            if k4mode=='mean': k4=np.mean(g4)*s4
            else:
                Xb=np.column_stack([np.ones(n),t*t]); cf=np.linalg.lstsq(Xb,g4,rcond=None)[0]; k4=(Xb@cf)*s4
            if k21_mc and rank>0:
                K21s=(zc*zc).T@zc/N_MC; np.fill_diagonal(K21s,0); R=K21s-K21m; np.fill_diagonal(R,0)
                Uu,sv,Vt=np.linalg.svd(R,full_matrices=False); K21=K21m+(Uu[:,:rank]*sv[:rank])@Vt[:rank]
            elif k21_mu:
                # inherited two-point cumulants along the known mean direction: K21_ab ~ u_a muhat_b, u = Cov(zt^2, A), A = muhat . zt
                muh=mu/np.linalg.norm(mu); A=zc@muh; u=((zc*zc).T@A)/N_MC; um=K21m@muh; K21=K21m+np.outer(u-um,muh)
                if k3_mu:
                    k3mu=(u-um)*muh*3.0   # diagonal of the rank-one slice: T_aaa = u_a mu_a (three permutations collapse) -- scale fitted below
                    Xb=np.column_stack([k3mu]); r=k3s-k3; cf=np.linalg.lstsq(Xb,r,rcond=None)[0]; k3=k3+cf[0]*k3mu
        t2=t*t; he2=t2-1; he4=(t2-6)*t2+3
        dE1=-k3*t*ph/(6*sz2)+(k3*k3)*he4*ph/(72*sz2*sz2*sz)+k4*he2*ph/(24*sz2*sz)
        dE2=k3*ph/(3*sz)+(k3*k3)*(3*t-t*t2)*ph/(36*sz2*sz2)-k4*t*ph/(12*sz2)
        m_next=Psi+dE1
        so=np.outer(sz,sz); rho=Sz/so; ds=T.hermite_d(t,ph,Ph,KH); rk=rho.copy(); C=np.zeros((n,n)); fact=1.0
        for k in range(1,KH+1):
            fact*=k; C+=rk*np.outer(ds[k-1],ds[k-1])/fact
            if k<KH: rk*=rho
        C*=so
        if k21_model or k21_mc:
            F1=Ph; F2=ph/sz; C+=0.5*(K21*np.outer(F2,F1)+K21.T*np.outer(F1,F2))
        var=(m2+dE2)-m_next*m_next if var_corr else m2-Psi*Psi
        C[np.diag_indices(n)]=var
        if K>0:
            M=Sz*Ph[:,None]; M[np.diag_indices(n)]=0; sources[l]=(M,ph/sz,2*Psi*(1-Ph),T.relu_k3(mu,sz,sz2,Ph,ph,Psi,m2),Ph)
        m=m_next; Sig=C; outs.append(m)
    return np.array(outs)
if __name__=='__main__':
    seeds,Y=T.load_nets(); nets=[int(a) for a in sys.argv[1].split(',')] if len(sys.argv)>1 else [0]
    configs=[dict(),dict(k21_mc=False),dict(k21_model=False),dict(rank=4),dict(rank=48),dict(K=1),dict(K=5),dict(N_MC=8000),dict(shrink_k3=0.0),dict(k4mode='fit')]
    for k in nets:
        W=T.regen(seeds[k]); Yk=Y[k]
        for cfg in configs:
            t0=time.time(); P=predict(W,**cfg); mse=np.mean((P-Yk)**2,axis=1)
            print(f"net {k} {cfg}: final={mse[-1]:.2e} [{time.time()-t0:.0f}s] ", " ".join(f"{x:.1e}" for x in mse[1:]), flush=True)
