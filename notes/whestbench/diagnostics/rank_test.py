# Is the far part of the cumulant path sum low rank?
# The contribution of source layer s to kappa_3(z_l) is a cubic form in the transport U^{(s->l)}
# (and N = M^T U).  Transport is right-multiplication by the same one-layer linear-response map for
# both, so a rank-r factorisation U = L R is preserved: only R (r x n) needs advancing.
# Test: keep the newest K0 sources exact, truncate older transports to rank r, compare with exact.
import sys, time, numpy as np
from scipy.special import ndtr
import twin2 as T
n=1024; L=16; INV=T.INV; pdf=T.pdf
def predict(W, N_MC=2400, K=8, K0=2, rank=0, T_ACT=2.75, KH=8, seed=0):
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
        k3=np.zeros(n); s_lo=max(l-K,0); U=Wl.astype(np.float64); age=0
        for src in range(l-1,s_lo-1,-1):
            Ms,ps,qs,k3hs,d1s=sources[src]
            if rank>0 and age>=K0:                      # finite-resolution transport of the far part
                Uu,sv,Vt=np.linalg.svd(U,full_matrices=False); U=(Uu[:,:rank]*sv[:rank])@Vt[:rank]
            Nt=Ms.T@U; UN=U*Nt
            k3+=((ps[:,None]*3)*UN*Nt+(qs[:,None]*3)*UN*U+k3hs[:,None]*(U*U*U)).sum(0)
            if src>s_lo: U=W[src]@(d1s[:,None]*U); age+=1
        z=Zs[l].astype(np.float64); a=z.mean(0); zc=z-a; v=(zc*zc).mean(0)
        k3s=(zc**3).mean(0); k4s=(zc**4).mean(0)-3*v*v
        r=(k3s-k3)/s3; Xb=np.column_stack([t,np.ones(n)]); cf=np.linalg.lstsq(Xb,r,rcond=None)[0]; k3=k3+(Xb@cf)*s3
        resid=k3s-k3; noise=6*np.mean(s3*s3)/N_MC; sig=max(np.mean(resid**2)-noise,0); k3=k3+(sig/(sig+noise))*resid
        k4=np.mean(k4s/s4)*s4
        muh=mu/np.linalg.norm(mu); Ap=zc@muh; u=(zc*zc).T@Ap/N_MC
        K21=np.zeros((n,n))  # exact source two-point slice
        t2=t*t; he2=t2-1; he4=(t2-6)*t2+3
        dE1=-k3*t*ph/(6*sz2)+(k3*k3)*he4*ph/(72*s4*sz)+k4*he2*ph/(24*s3)
        dE2=k3*ph/(3*sz)+(k3*k3)*(3*t-t*t2)*ph/(36*s4)-k4*t*ph/(12*sz2)
        m_next=Psi+dE1
        so=np.outer(sz,sz); rho=Sz/so; ds=T.hermite_d(t,ph,Ph,KH); rk=rho.copy(); C=np.zeros((n,n)); fact=1.0
        for k in range(1,KH+1):
            fact*=k; C+=rk*np.outer(ds[k-1],ds[k-1])/fact
            if k<KH: rk*=rho
        C*=so
        F2=ph/sz; K21=np.outer(u,muh); np.fill_diagonal(K21,0)
        C+=0.5*(K21*np.outer(F2,Ph)+K21.T*np.outer(Ph,F2))
        C[np.diag_indices(n)]=(m2+dE2)-m_next*m_next
        Msrc=Sz*Ph[:,None]; Msrc[np.diag_indices(n)]=0
        m3=mu*m2+2*sz2*Psi; sources[l]=(Msrc,ph/sz,2*Psi*(1-Ph),m3-3*m2*Psi+2*Psi**3,Ph)
        m=m_next; Sig=C; outs.append(m)
    return np.array(outs)
if __name__=='__main__':
    seeds,Y=T.load_nets()
    for k in [0,1]:
        W=T.regen(seeds[k]); Yk=Y[k]
        for cfg in [dict(K=1),dict(K=8,rank=0),dict(K=8,K0=2,rank=128),dict(K=8,K0=2,rank=64),dict(K=8,K0=2,rank=32),dict(K=8,K0=1,rank=64),dict(K=12,K0=2,rank=64),dict(K=8,K0=3,rank=48)]:
            t0=time.time(); P=predict(W,**cfg); mse=np.mean((P-Yk)**2,axis=1)
            print(f"net {k} {cfg}: final={mse[-1]:.3e} [{time.time()-t0:.0f}s]", flush=True)
