# v2 numpy twin: Gaussian-closure propagation + exact recent k3 sources + on-the-fly MC calibration of the
# aligned (few-parameter) closure residuals per layer. Evaluated on shard-0 networks with 1e9 truth.
import sys, time, numpy as np, pyarrow.parquet as pq, pyarrow as pa
from scipy.special import ndtr
n=1024; L=16; SQ2PI=np.sqrt(2*np.pi); INV=1/SQ2PI
pdf=lambda t: np.exp(-0.5*t*t)*INV
def load_nets():
    pf=pq.ParquetFile('../data/mini-00000-of-00007.parquet'); tb=pf.read_row_group(0,columns=['mlp_seed','all_layer_means'])
    seeds=tb.column('mlp_seed').to_pylist(); a=tb.column('all_layer_means').combine_chunks()
    while pa.types.is_list(a.type) or pa.types.is_fixed_size_list(a.type) or pa.types.is_large_list(a.type): a=a.flatten()
    Y=a.to_numpy(zero_copy_only=False).astype(np.float64).reshape(-1,L,n); return seeds,Y
def regen(seed):
    ss=np.random.SeedSequence(int(seed)).spawn(3); rng=np.random.default_rng(ss[0]); sc=float(np.sqrt(2.0/n))
    return [(rng.standard_normal((n,n))*sc).astype(np.float32) for _ in range(L)]
def hermite_d(t,ph,Ph,kmax):
    x=-t; ds=[Ph,ph]; hp,h=np.ones_like(x),x
    for k in range(3,kmax+1):
        ds.append(ph*h); hp,h=h,x*h-(k-2)*hp
    return ds[:kmax]
def relu_k3(mu,s,s2,Ph,ph,m1,m2):
    m3=mu*m2+2*s2*m1; return m3-3*m2*m1+2*m1**3
def predict(W, N_MC=4000, K=3, T_ACT=3.0, KH=8, cal_mean=True, cal_var=True, seed=0, use_k3sq=True, var_k3=False, basis='e3'):
    rng=np.random.default_rng(seed)
    # MC samples (float32), per-layer pre-activation samples kept
    if N_MC>0:
        X=rng.standard_normal((N_MC,n)).astype(np.float32); Zs=[]; h=X
        for l in range(L):
            z=h@W[l]; Zs.append(z); h=np.maximum(z,0)
    outs=[]; sources={}
    W0=W[0]; G=W0.T@W0; s2=np.diag(G); s=np.sqrt(s2); m=s*INV; so=np.outer(s,s)
    rho=np.clip(G/so,-1,1); C=so*((np.sqrt(1-rho*rho)+rho*(np.pi-np.arccos(rho)))/(2*np.pi))-np.outer(m,m)
    C[np.diag_indices(n)]=s2*(0.5-1/(2*np.pi))
    M=G*0.5; M[np.diag_indices(n)]=0
    sources[0]=(M,INV/s,m,s2*s*(2*INV-1.5*INV+2*INV**3),np.full(n,0.5,dtype=np.float32))
    outs.append(m); Sig=C
    for l in range(1,L):
        Wl=W[l]; mu=m@Wl; Sz=(Wl.T@Sig)@Wl; sz2=np.diag(Sz).copy()
        if N_MC>0 and cal_var:
            # calibrate propagated pre-activation variances: regress sample var(z)/sz2 - 1 on (1,t,t^2)
            zl=Zs[l].astype(np.float64); vs=zl.var(0); t0=mu/np.sqrt(sz2)
            Xb=np.column_stack([np.ones(n),t0,t0*t0]); cf=np.linalg.lstsq(Xb,vs/sz2-1,rcond=None)[0]
            fac=1+Xb@cf; sz2=sz2*fac; Sz=Sz*np.sqrt(np.outer(fac,fac))
        sz=np.sqrt(sz2); t=mu/sz; ph=pdf(t); Ph=ndtr(t); Psi=mu*Ph+sz*ph; m2=(mu*mu+sz2)*Ph+mu*sz*ph
        k3=np.zeros(n)
        if K>0:
            act=np.nonzero(np.abs(t)<T_ACT)[0]; s_lo=max(l-K,0); U=Wl[:,act].astype(np.float64); k3a=np.zeros(act.size)
            for src in range(l-1,s_lo-1,-1):
                Ms,ps,qs,k3hs,d1s=sources[src]; Nt=Ms.T@U; UN=U*Nt
                k3a+=((ps[:,None]*3)*UN*Nt+(qs[:,None]*3)*UN*U+k3hs[:,None]*(U*U*U)).sum(0)
                if src>s_lo: U=W[src]@(d1s[:,None]*U)
            k3[act]=k3a
        tph=t*ph; corr=-k3*tph/(6*sz2)
        if use_k3sq: t2=t*t; corr=corr+(k3*k3)*((t2-6)*t2+3)*ph/(72*sz2*sz2*sz)
        m_next=Psi+corr
        if N_MC>0 and cal_mean:
            # residual target: sample mean of relu(z) minus Gaussian closure at the SAMPLE moments (quadratic control)
            zl=Zs[l].astype(np.float64); mz=zl.mean(0); vz=zl.var(0); szh=np.sqrt(vz); th=mz/szh
            r=np.maximum(zl,0).mean(0)-(mz*ndtr(th)+szh*pdf(th))-corr   # remove what the exact channel already explains
            if basis=='e3': Xb=np.column_stack([sz*ph*t, sz*ph, sz*ph*(t*t-1)])
            else: Xb=np.column_stack([sz*ph*t, sz*ph, sz*ph*(t*t-1), sz*ph*t*(t*t-3)])
            cf=np.linalg.lstsq(Xb,r,rcond=None)[0]; m_next=m_next+Xb@cf
        so=np.outer(sz,sz); rho=Sz/so; ds=hermite_d(t,ph,Ph,KH); rk=rho.copy(); C=np.zeros((n,n)); fact=1.0
        for k in range(1,KH+1):
            fact*=k; dk=ds[k-1]; C+=rk*np.outer(dk,dk)/fact
            if k<KH: rk*=rho
        C*=so; var=m2-Psi*Psi
        if var_k3: var=var+k3*ph/(3*sz)-2*Psi*(m_next-Psi)
        C[np.diag_indices(n)]=var
        if K>0:
            M=Sz*Ph[:,None]; M[np.diag_indices(n)]=0; sources[l]=(M,ph/sz,2*Psi*(1-Ph),relu_k3(mu,sz,sz2,Ph,ph,Psi,m2),Ph)
        m=m_next; Sig=C; outs.append(m)
    return np.array(outs)
if __name__=='__main__':
    seeds,Y=load_nets(); nets=[int(a) for a in sys.argv[1].split(',')] if len(sys.argv)>1 else [0]
    configs=[dict(N_MC=0,K=0),dict(N_MC=0,K=3),dict(N_MC=4000,K=0,cal_mean=False),dict(N_MC=4000,K=0,cal_var=False),dict(N_MC=4000,K=0),dict(N_MC=4000,K=3),dict(N_MC=4000,K=3,basis='e4'),dict(N_MC=8000,K=3)]
    for k in nets:
        W=regen(seeds[k]); Yk=Y[k]
        for cfg in configs:
            t0=time.time(); P=predict(W,**cfg); mse=np.mean((P-Yk)**2,axis=1)
            print(f"net {k} {cfg}: final={mse[-1]:.2e} [{time.time()-t0:.0f}s] ", " ".join(f"{x:.1e}" for x in mse[1:]), flush=True)
