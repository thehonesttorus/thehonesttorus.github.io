# v3 twin: Gaussian-closure propagation + exact recent k3 sources + MC-fitted aligned cumulant models
#   k3_j = exact_recent_j + (a t_j + b) s_j^3   [a,b regressed over neurons on sample k3 - exact_recent]
#   k4_j = (c + d t_j^2 + e t_j) s_j^4            [regressed over neurons on sample k4]
# readout E4 (k3, k3^2, k4 terms); variance closure with k3, k3^2, k4 (diagonal); Gaussian off-diagonal.
import sys, time, numpy as np
from scipy.special import ndtr
import twin2 as T
n=1024; L=16; INV=T.INV; pdf=T.pdf
def cums_from_samples(z):
    z=z.astype(np.float64); a=z.mean(0); zc=z-a; v=(zc*zc).mean(0); k3=(zc**3).mean(0); k4=(zc**4).mean(0)-3*v*v
    return a,v,k3,k4
def predict(W, N_MC=4000, K=3, T_ACT=3.0, KH=8, seed=0, fit_k3=True, fit_k4=True, var_corr=True, shrink=0.0, k4basis=1, k3basis=2, oracle_cov=None, oracle_mean=None):
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
        Psi=mu*Ph+sz*ph; m2=(mu*mu+sz2)*Ph+mu*sz*ph
        k3=np.zeros(n)
        if K>0:
            act=np.nonzero(np.abs(t)<T_ACT)[0]; s_lo=max(l-K,0); U=Wl[:,act].astype(np.float64); k3a=np.zeros(act.size)
            for src in range(l-1,s_lo-1,-1):
                Ms,ps,qs,k3hs,d1s=sources[src]; Nt=Ms.T@U; UN=U*Nt
                k3a+=((ps[:,None]*3)*UN*Nt+(qs[:,None]*3)*UN*U+k3hs[:,None]*(U*U*U)).sum(0)
                if src>s_lo: U=W[src]@(d1s[:,None]*U)
            k3[act]=k3a
        k4=np.zeros(n)
        if Zs is not None:
            a_,v_,k3s,k4s=cums_from_samples(Zs[l]); s3=sz2*sz; s4=sz2*sz2
            if fit_k3:
                r=(k3s-k3)/s3
                Xb=np.column_stack([t,np.ones(n)]) if k3basis==2 else np.column_stack([t,np.ones(n),t*t])
                cf=np.linalg.lstsq(Xb,r,rcond=None)[0]; k3=k3+(Xb@cf)*s3
                if shrink>0: k3=k3+shrink*(k3s-k3)
            if fit_k4:
                r4=k4s/s4
                Xb=np.column_stack([np.ones(n)]) if k4basis==1 else np.column_stack([np.ones(n),t*t,t])
                cf=np.linalg.lstsq(Xb,r4,rcond=None)[0]; k4=(Xb@cf)*s4
        t2=t*t; he2=t2-1; he4=(t2-6)*t2+3
        dE1=-k3*t*ph/(6*sz2)+(k3*k3)*he4*ph/(72*sz2*sz2*sz)+k4*he2*ph/(24*sz2*sz)
        dE2=k3*ph/(3*sz)+(k3*k3)*(3*t-t*t2)*ph/(36*sz2*sz2)-k4*t*ph/(12*sz2)
        m_next=Psi+dE1
        so=np.outer(sz,sz); rho=Sz/so; ds=T.hermite_d(t,ph,Ph,KH); rk=rho.copy(); C=np.zeros((n,n)); fact=1.0
        for k in range(1,KH+1):
            fact*=k; C+=rk*np.outer(ds[k-1],ds[k-1])/fact
            if k<KH: rk*=rho
        C*=so; var=(m2+dE2)-m_next*m_next if var_corr else m2-Psi*Psi
        C[np.diag_indices(n)]=var
        if K>0:
            M=Sz*Ph[:,None]; M[np.diag_indices(n)]=0; sources[l]=(M,ph/sz,2*Psi*(1-Ph),T.relu_k3(mu,sz,sz2,Ph,ph,Psi,m2),Ph)
        m=m_next; Sig=C; outs.append(m)
    return np.array(outs)
if __name__=='__main__':
    seeds,Y=T.load_nets(); nets=[int(a) for a in sys.argv[1].split(',')] if len(sys.argv)>1 else [0]
    configs=[dict(N_MC=4000,K=3),dict(N_MC=4000,K=3,var_corr=False),dict(N_MC=4000,K=3,fit_k4=False),dict(N_MC=4000,K=1),dict(N_MC=4000,K=5),dict(N_MC=4000,K=3,k4basis=3),dict(N_MC=4000,K=3,k3basis=3),dict(N_MC=4000,K=3,shrink=0.3),dict(N_MC=8000,K=3),dict(N_MC=2000,K=3)]
    for k in nets:
        W=T.regen(seeds[k]); Yk=Y[k]
        for cfg in configs:
            t0=time.time(); P=predict(W,**cfg); mse=np.mean((P-Yk)**2,axis=1)
            print(f"net {k} {cfg}: final={mse[-1]:.2e} [{time.time()-t0:.0f}s] ", " ".join(f"{x:.1e}" for x in mse[1:]), flush=True)
