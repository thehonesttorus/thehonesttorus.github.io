# How much depth do we need, GIVEN that a fitted kernel absorbs the aggregate?
# Caches, per network and layer: t, sigma^3, every partial source sum S_j (j=1..l), and a
# high-accuracy MC third cumulant.  Everything downstream is then a linear solve.
import numpy as np, sys, os
from scipy.special import ndtr
import twin2 as T
n=1024; L=16; INV=T.INV; pdf=T.pdf
seeds,Yall=T.load_nets()
SCHED=np.array([[0.00094,0.00715],[0.00783,0.01340],[0.01134,0.01947],[0.01406,0.02302],[0.01658,0.02679],
 [0.01853,0.03065],[0.02040,0.03478],[0.02217,0.03824],[0.02331,0.04162],[0.02464,0.04615],
 [0.02555,0.04940],[0.02684,0.05096],[0.02781,0.05280],[0.02838,0.05462],[0.02940,0.05873]])
C2P=np.array([0.00396,0.07819,0.14538,0.22329,0.30105,0.38174,0.45643,0.55252,0.63162,0.70963,
              0.80605,0.92641,1.00444,1.09326,1.18779])
NB=int(sys.argv[2]) if len(sys.argv)>2 else 60000
NETS=[int(x) for x in sys.argv[1].split(',')]

def run(k, N_MC=2400, seed=0):
    W=T.regen(seeds[k]); rng=np.random.default_rng(seed)
    X=rng.standard_normal((N_MC,n)).astype(np.float32); Zs=[]; h=X
    for l in range(L): z=h@W[l]; Zs.append(z); h=np.maximum(z,0)
    rb=np.random.default_rng(9000+k); S1=np.zeros((L,n)); S2=np.zeros((L,n)); S3=np.zeros((L,n)); ch=5000; done=0
    while done<NB:
        b=min(ch,NB-done); Xb=rb.standard_normal((b,n)).astype(np.float32); hb=Xb
        for l in range(L):
            z=(hb@W[l]).astype(np.float64); S1[l]+=z.sum(0); S2[l]+=(z*z).sum(0); S3[l]+=(z**3).sum(0)
            hb=np.maximum(z,0).astype(np.float32)
        done+=b
    M1=S1/NB; M2=S2/NB; M3=S3/NB; K3T=M3-3*M2*M1+2*M1**3
    out={}; sources={}
    W0=W[0]; G=W0.T@W0; s2=np.diag(G); s=np.sqrt(s2); m=s*INV; so=np.outer(s,s)
    rho=np.clip(G/so,-1,1); C=so*((np.sqrt(1-rho*rho)+rho*(np.pi-np.arccos(rho)))/(2*np.pi))-np.outer(m,m)
    C[np.diag_indices(n)]=s2*(0.5-1/(2*np.pi)); M=G*0.5; M[np.diag_indices(n)]=0
    sources[0]=(M,INV/s,m,s2*s*(2*INV-1.5*INV+2*INV**3),np.full(n,0.5)); Sig=C
    for l in range(1,L):
        Wl=W[l]; mu=m@Wl; Sz=(Wl.T@Sig)@Wl; sz2=np.diag(Sz).copy(); sz=np.sqrt(np.maximum(sz2,1e-30)); t=mu/sz
        ph=pdf(t); Ph=ndtr(t); Psi=mu*Ph+sz*ph; m2=(mu*mu+sz2)*Ph+mu*sz*ph; s3=sz2*sz; s4=sz2*sz2
        U=Wl.astype(np.float64); part=[]; run_=np.zeros(n)
        for srcl in range(l-1,-1,-1):
            Ms,ps,qs,k3hs,d1s=sources[srcl]; Nt=Ms.T@U; UN=U*Nt
            run_=run_+((ps[:,None]*3)*UN*Nt+(qs[:,None]*3)*UN*U+k3hs[:,None]*(U*U*U)).sum(0)
            part.append(run_.copy())
            if srcl>0: U=W[srcl]@(d1s[:,None]*U)
        out[l]=dict(t=t.astype(np.float32),s3=s3.astype(np.float32),
                    part=np.array(part,dtype=np.float32),truth=K3T[l].astype(np.float32),
                    Y=Yall[k][l].astype(np.float32),sz=sz.astype(np.float32))
        k3_1=part[0]
        Ud=Wl.astype(np.float64); Ms,ps,qs,k3hs,d1s=sources[l-1]; Ntd=Ms.T@Ud; UNd=Ud*Ntd
        X1=(UNd*2)*ps[:,None]+(Ud*Ud)*qs[:,None]
        X2=(Ntd*Ntd)*ps[:,None]+(UNd*2)*qs[:,None]+(Ud*Ud)*k3hs[:,None]
        K21=X1.T@Ntd+X2.T@Ud; np.fill_diagonal(K21,0.0)
        z=Zs[l].astype(np.float64); a=z.mean(0); zc=z-a; k3s=(zc**3).mean(0)
        k3=k3_1+SCHED[l-1,0]*t*s3; g4=SCHED[l-1,1]
        resid=k3s-k3; noise=6*np.mean(s3*s3)/N_MC; sig=max(np.mean(resid**2)-noise,0)
        k3=k3+(sig/(sig+noise))*resid
        muh=mu/np.linalg.norm(mu); Ap=zc@muh
        u=0.75*(C2P[l-1]*s3)+0.25*((zc*zc).T@Ap/N_MC-K21@muh)
        K21=K21+np.outer(u,muh); np.fill_diagonal(K21,0.0)
        k4v=g4*s4; t2=t*t; he2=t2-1; he4=(t2-6)*t2+3
        dE1=-k3*t*ph/(6*sz2)+(k3*k3)*he4*ph/(72*s4*sz)+k4v*he2*ph/(24*s3)
        dE2=k3*ph/(3*sz)+(k3*k3)*(3*t-t*t2)*ph/(36*s4)-k4v*t*ph/(12*sz2)
        m_next=Psi+dE1
        so=np.outer(sz,sz); rhoc=Sz/so; ds=T.hermite_d(t,ph,Ph,8); rk=rhoc.copy(); Cn=np.zeros((n,n)); fact=1.0
        for kk in range(1,9):
            fact*=kk; Cn+=rk*np.outer(ds[kk-1],ds[kk-1])/fact
            if kk<8: rk*=rhoc
        Cn*=so; F2=ph/sz
        Cn=Cn+0.5*(K21*np.outer(F2,Ph)+K21.T*np.outer(Ph,F2))
        Cn[np.diag_indices(n)]=(m2+dE2)-m_next*m_next
        Msrc=Sz*Ph[:,None]; np.fill_diagonal(Msrc,0.0); m3=mu*m2+2*sz2*Psi
        sources[l]=(Msrc,ph/sz,2*Psi*(1-Ph),m3-3*m2*Psi+2*Psi**3,Ph)
        m=m_next; Sig=Cn
    return out

for k in NETS:
    fn=f"depthcache_{k}.npz"
    if os.path.exists(fn): print("skip",k,flush=True); continue
    o=run(k); d={}
    for l in range(1,L):
        for key,v in o[l].items(): d[f"{l}_{key}"]=v
    np.savez_compressed(fn,**d); print("cached",k,flush=True)
