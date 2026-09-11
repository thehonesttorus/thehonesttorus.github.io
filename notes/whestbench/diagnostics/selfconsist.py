# A free self-consistent source correction.  The source parameters handed to the next layer are
# computed from a GAUSSIAN z.  But we have a model of kappa_3(z).  The exact Hermite integral
#     int_{-t}^{inf} (x+t)^m phi(x) He_k(x) dx = m! phi(t) He_{k-m-1}(-t)   (m<k),   J_3^3 = 6 Phi(t)
# gives the first-order response of the ReLU moments to the input's third cumulant:
#     dE[h]   = -k3 t phi / (6 s^2)      (this is already the Edgeworth readout)
#     dE[h^2] =  k3 phi / (3 s)          (this is already the variance correction)
#     dE[h^3] =  k3 Phi
# hence   d kappa_3(h) = k3 * [ Phi - phi Psi / s + m2 t phi/(2 s^2) - Psi^2 t phi / s^2 ].
# Adding that to the source's diagonal entry makes the depth-1 channel carry inherited
# non-Gaussianity recursively, at O(n) cost -- the "induce the higher correction from the lower one"
# that a Leibniz-coherent closure would do, instead of fitting it separately.
import numpy as np, sys
from scipy.special import ndtr
import ablate as AB
import twin2 as T
n=1024; L=16; INV=T.INV; pdf=T.pdf
seeds,Yall=AB.seeds,AB.Yall
A_S,G_S,C_S=AB.A_S,AB.G_S,AB.C_S

def predict(k,cfg,seed=0):
    W=AB.getW(k); N=cfg['N']; dep=cfg['depth']; k21x=cfg['k21']
    al,ga,de,be=cfg['al'],cfg['ga'],cfg['de'],cfg.get('be',0.0)
    bm=cfg.get('bm',0.0)
    if N>0:
        rng=np.random.default_rng(seed); X=rng.standard_normal((N,n)).astype(np.float32); Zs=[]; h=X
        for l in range(L): z=h@W[l]; Zs.append(z); h=np.maximum(z,0)
    sources={}
    W0=W[0]; G=W0.T@W0; s2=np.diag(G); s=np.sqrt(s2); m=s*INV; so=np.outer(s,s)
    rho=np.clip(G/so,-1,1); C=so*((np.sqrt(1-rho*rho)+rho*(np.pi-np.arccos(rho)))/(2*np.pi))-np.outer(m,m)
    C[np.diag_indices(n)]=s2*(0.5-1/(2*np.pi)); C=C.astype(np.float32)
    M=(G*0.5).astype(np.float32); M[np.diag_indices(n)]=0
    sources[0]=(M,(INV/s).astype(np.float32),m.astype(np.float32),
                (s2*s*(2*INV-1.5*INV+2*INV**3)).astype(np.float32),np.full(n,0.5,np.float32))
    Sig=C; m=m.astype(np.float32)
    for l in range(1,L):
        Wl=W[l]; mu=(m@Wl).astype(np.float64); Sz=((Wl.T@Sig)@Wl).astype(np.float64)
        sz2=np.diag(Sz).copy(); sz=np.sqrt(np.maximum(sz2,1e-30)); t=mu/sz
        ph=pdf(t); Ph=ndtr(t); Psi=mu*Ph+sz*ph; m2=(mu*mu+sz2)*Ph+mu*sz*ph; s3=sz2*sz; s4=sz2*sz2
        k3=np.zeros(n,np.float32); K21=np.zeros((n,n),np.float32); d=dep(l)
        if d>0:
            U=Wl; lo=max(l-d,0)
            for srcl in range(l-1,lo-1,-1):
                if srcl not in sources: break
                Ms,ps,qs,k3hs,d1s=sources[srcl]; Nt=Ms.T@U; UN=U*Nt
                k3=k3+((ps[:,None]*np.float32(3))*UN*Nt+(qs[:,None]*np.float32(3))*UN*U
                       +k3hs[:,None]*(U*U*U)).sum(0)
                if k21x:
                    X1=(UN*2)*ps[:,None]+(U*U)*qs[:,None]
                    X2=(Nt*Nt)*ps[:,None]+(UN*2)*qs[:,None]+(U*U)*k3hs[:,None]
                    Kc=X1.T@Nt+X2.T@U; np.fill_diagonal(Kc,0.0); K21=K21+Kc
                if srcl>lo: U=W[srcl]@(d1s[:,None]*U)
        k3=k3.astype(np.float64)+al*A_S[l-1]*t*s3; g4=de*G_S[l-1]
        muh=mu/np.linalg.norm(mu)
        if N>0:
            z=Zs[l]; a=z.mean(0); zc=z-a; k3s=(zc**3).mean(0).astype(np.float64)
            resid=k3s-k3; noise=6*np.mean(s3*s3)/N; sg=max(np.mean(resid**2)-noise,0)
            k3=k3+(sg/(sg+noise))*resid
            Ap=zc@muh.astype(np.float32); u_mc=((zc*zc).T@Ap).astype(np.float64)/N-K21@muh
            u=0.75*(ga*C_S[l-1]*s3)+0.25*u_mc
        else:
            u=ga*C_S[l-1]*s3
        K21=K21.astype(np.float64)+np.outer(u,muh); np.fill_diagonal(K21,0.0)
        k4v=g4*s4; t2=t*t; he2=t2-1; he4=(t2-6)*t2+3
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
        # ---- source state, with the free self-consistent corrections ----
        Msrc=Sz*Ph[:,None]; np.fill_diagonal(Msrc,0.0)
        m3=mu*m2+2*sz2*Psi
        k3h=m3-3*m2*Psi+2*Psi**3
        if be!=0.0:
            k3h=k3h+be*k3*(Ph-ph*Psi/sz+m2*t*ph/(2*sz2)-Psi*Psi*t*ph/sz2)
        if bm!=0.0:
            # two-point analogue: the source matrix M_{ca}=Sz_{ca}Phi(t_c) picks up the tracked
            # two-point cumulant kappa(z_c,z_c,z_a) with the same first-order response coefficient
            Msrc=Msrc+bm*(K21*(ph/sz)[:,None])
            np.fill_diagonal(Msrc,0.0)
        sources[l]=(Msrc.astype(np.float32),(ph/sz).astype(np.float32),
                    (2*Psi*(1-Ph)).astype(np.float32),k3h.astype(np.float32),Ph.astype(np.float32))
        m=m_next.astype(np.float32); Sig=Cn.astype(np.float32)
    return np.mean((m-Yall[k][-1])**2)

FIT=[0,3,6]; HOLD=list(range(8,16))
def sc(cfg,nets): return float(np.mean([predict(k,cfg) for k in nets]))
base=dict(AB.CFG['A deployed']); base['be']=0.0; base['bm']=0.0
print("Self-consistent source corrections on the DEPLOYED configuration (cost unchanged, 198 n^3).\n")
print(f"{'beta_k3h':>9} {'beta_M':>7} {'alpha':>6}   {'fit MSE':>10} {'held MSE':>10}   note")
for be,bm in [(0.0,0.0),(0.5,0.0),(1.0,0.0),(1.5,0.0),(2.0,0.0),(0.0,1.0),(1.0,1.0),(1.0,0.5),(1.0,-1.0)]:
    c=dict(base); c['be']=be; c['bm']=bm
    best=dict(c); bs=sc(best,FIT)
    for g in (0.4,0.6,0.8,1.2):
        cc=dict(c); cc['al']=g; s=sc(cc,FIT)
        if s<bs: bs=s; best=cc
    hs=sc(best,HOLD)
    tag="theory value" if (be==1.0 and bm in (0.0,1.0)) else ""
    print(f"{be:9.2f} {bm:7.2f} {best['al']:6.2f}   {bs:10.3e} {hs:10.3e}   {tag}", flush=True)
