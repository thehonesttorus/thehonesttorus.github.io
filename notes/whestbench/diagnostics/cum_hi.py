# The k3 oracle was contaminated.  lambda_3 = kappa_3/sigma^3 is only about 0.03 for these networks,
# while the Monte-Carlo error on it is sqrt(6/N) = 0.0071 at N=120000 -- 24% noise.  Its induced
# readout MSE is  mean[ (sigma t phi/6)^2 * 6/N + (sigma (t^2-1) phi/24)^2 * 24/N ],  which is of the
# same order as the residual we were trying to measure.  Redo at N = 1e6 and subtract the residual
# noise analytically.
import numpy as np, sys, os
from scipy.special import ndtr
import ablate as AB, twin2 as T
n=1024; L=16; INV=T.INV; pdf=T.pdf
seeds,Yall=AB.seeds,AB.Yall
A_S,G_S,C_S=AB.A_S,AB.G_S,AB.C_S
NB=int(sys.argv[1]) if len(sys.argv)>1 else 1000000
NETS=[int(x) for x in sys.argv[2].split(',')] if len(sys.argv)>2 else [0,1]

def moments(k,NB):
    fn=f"mom_{k}_{NB}.npz"
    if os.path.exists(fn):
        d=np.load(fn); return d['lam3'],d['lam4'],d['sd']
    W=AB.getW(k); rb=np.random.default_rng(555+k); ch=4000; done=0
    S=[np.zeros((L,n)) for _ in range(4)]
    while done<NB:
        b=min(ch,NB-done); Xb=rb.standard_normal((b,n)).astype(np.float32); hb=Xb
        for l in range(L):
            z=(hb@W[l]).astype(np.float64); z2=z*z
            S[0][l]+=z.sum(0); S[1][l]+=z2.sum(0); S[2][l]+=(z2*z).sum(0); S[3][l]+=(z2*z2).sum(0)
            hb=np.maximum(z,0).astype(np.float32)
        done+=b
        if done % 100000 < ch: print("  net",k,done,flush=True)
    m1,m2,m3,m4=[x/NB for x in S]
    v=m2-m1*m1; sd=np.sqrt(v)
    k3=m3-3*m2*m1+2*m1**3
    k4=m4-4*m3*m1-3*m2*m2+12*m2*m1*m1-6*m1**4
    lam3=k3/sd**3; lam4=k4/sd**4
    np.savez_compressed(fn,lam3=lam3,lam4=lam4,sd=sd); return lam3,lam4,sd

def run(k,lam,exact_mean,oracle,N=2400,seed=0):
    W=AB.getW(k); Y=Yall[k]
    rng=np.random.default_rng(seed); X=rng.standard_normal((N,n)).astype(np.float32); Zs=[]; h=X
    for l in range(L): z=h@W[l]; Zs.append(z); h=np.maximum(z,0)
    sources={}; noise_mse=[]
    W0=W[0]; G=W0.T@W0; s2=np.diag(G); s=np.sqrt(s2); m=s*INV; so=np.outer(s,s)
    rho=np.clip(G/so,-1,1); C=so*((np.sqrt(1-rho*rho)+rho*(np.pi-np.arccos(rho)))/(2*np.pi))-np.outer(m,m)
    C[np.diag_indices(n)]=s2*(0.5-1/(2*np.pi)); C=C.astype(np.float32)
    M=(G*0.5).astype(np.float32); M[np.diag_indices(n)]=0
    sources[0]=(M,(INV/s).astype(np.float32),m.astype(np.float32),
                (s2*s*(2*INV-1.5*INV+2*INV**3)).astype(np.float32),np.full(n,0.5,np.float32))
    Sig=C; m=((Y[0] if exact_mean else m)).astype(np.float32)
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
        rr=k3s-k3; nz=6*np.mean(s3*s3)/N; sg=max(np.mean(rr**2)-nz,0); k3=k3+(sg/(sg+nz))*rr
        k4v=g4*s4
        if 'k3' in oracle: k3=lam[0][l]*s3
        if 'k4' in oracle: k4v=lam[1][l]*s4
        noise_mse.append(np.mean((sz*t*ph/6)**2)*6/NB*(1 if 'k3' in oracle else 0)
                        +np.mean((sz*(t*t-1)*ph/24)**2)*24/NB*(1 if 'k4' in oracle else 0))
        muh=mu/np.linalg.norm(mu); Ap=zc@muh.astype(np.float32)
        u=0.75*(C_S[l-1]*s3)+0.25*(((zc*zc).T@Ap).astype(np.float64)/N-K21@muh)
        K21=K21.astype(np.float64)+np.outer(u,muh); np.fill_diagonal(K21,0.0)
        t2=t*t; he2=t2-1; he4=(t2-6)*t2+3
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
        m=((Y[l] if (exact_mean and l<L-1) else m_next)).astype(np.float32)
    return float(np.mean((m.astype(np.float64)-Y[-1])**2)), noise_mse[-1]

LAM={k:moments(k,NB) for k in NETS}
print(f"\nCumulant oracle at N = {NB}.  'raw' is measured; 'net' subtracts the analytic Monte-Carlo")
print("noise the substituted cumulants inject into the final readout.\n")
print(f"{'exact mean':>10} {'oracle':>10}   " + "  ".join(f"net {k}" for k in NETS) +
      f"   {'raw mean':>10} {'noise':>10} {'net':>10}")
for em in (False,True):
    for ora in ([], ['k3'], ['k3','k4']):
        res=[run(k,LAM[k],em,ora) for k in NETS]
        v=float(np.mean([r[0] for r in res])); nz=float(np.mean([r[1] for r in res]))
        if not em: nz=0.0  # noise at intermediate layers propagates; only report it for the isolated step
        print(f"{str(em):>10} {'+'.join(ora) or 'none':>10}   " + "  ".join(f"{r[0]:.3e}" for r in res) +
              f"   {v:.3e} {nz:10.2e} {max(v-nz,0):10.3e}", flush=True)
