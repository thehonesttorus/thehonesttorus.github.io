# What is the ceiling for improving the cubic channel?  Substitute the MEASURED standardised cumulants
# of the pre-activation into the deployed closure, layer by layer, and read the final MSE.
# Uses the Hermite-moment cache written by readout_wall.py:  Hm[l,k] = E[He_k((z-mu)/sigma)],
# so lambda_3 = Hm[l,3], lambda_4 = Hm[l,4], lambda_5 = Hm[l,5], lambda_6 = Hm[l,6] - 10*Hm[l,3]^2.
import numpy as np, sys, math, os
from scipy.special import ndtr
import twin2 as T
n=1024; L=16; INV=T.INV; pdf=T.pdf
seeds,Yall=T.load_nets()
SCHED=np.array([[0.00094,0.00715],[0.00783,0.01340],[0.01134,0.01947],[0.01406,0.02302],[0.01658,0.02679],
 [0.01853,0.03065],[0.02040,0.03478],[0.02217,0.03824],[0.02331,0.04162],[0.02464,0.04615],
 [0.02555,0.04940],[0.02684,0.05096],[0.02781,0.05280],[0.02838,0.05462],[0.02940,0.05873]])
C2P=np.array([0.00396,0.07819,0.14538,0.22329,0.30105,0.38174,0.45643,0.55252,0.63162,0.70963,
              0.80605,0.92641,1.00444,1.09326,1.18779])
NB=int(sys.argv[1]) if len(sys.argv)>1 else 120000
NETS=[int(x) for x in sys.argv[2].split(',')] if len(sys.argv)>2 else [0,1]

def he_all(x,K):
    H=[np.ones_like(x),x]
    for k in range(2,K+1): H.append(x*H[k-1]-(k-1)*H[k-2])
    return H

def predict(k, ora_k3=None, ora_k4=None, hi_order=False, N_MC=2400, seed=0, lam=None):
    """ora_k3/ora_k4: None, or a predicate on l saying 'use the measured standardised cumulant'."""
    W=T.regen(seeds[k]); rng=np.random.default_rng(seed)
    X=rng.standard_normal((N_MC,n)).astype(np.float32); Zs=[]; h=X
    for l in range(L): z=h@W[l]; Zs.append(z); h=np.maximum(z,0)
    sources={}; outs=[]
    W0=W[0]; G=W0.T@W0; s2=np.diag(G); s=np.sqrt(s2); m=s*INV; so=np.outer(s,s)
    rho=np.clip(G/so,-1,1); C=so*((np.sqrt(1-rho*rho)+rho*(np.pi-np.arccos(rho)))/(2*np.pi))-np.outer(m,m)
    C[np.diag_indices(n)]=s2*(0.5-1/(2*np.pi)); M=G*0.5; M[np.diag_indices(n)]=0
    sources[0]=(M,INV/s,m,s2*s*(2*INV-1.5*INV+2*INV**3),np.full(n,0.5)); Sig=C; outs.append(m)
    for l in range(1,L):
        Wl=W[l]; mu=m@Wl; Sz=(Wl.T@Sig)@Wl; sz2=np.diag(Sz).copy(); sz=np.sqrt(np.maximum(sz2,1e-30)); t=mu/sz
        ph=pdf(t); Ph=ndtr(t); Psi=mu*Ph+sz*ph; m2=(mu*mu+sz2)*Ph+mu*sz*ph; s3=sz2*sz; s4=sz2*sz2
        Ud=Wl.astype(np.float64); Ms,ps,qs,k3hs,d1s=sources[l-1]; Ntd=Ms.T@Ud; UNd=Ud*Ntd
        k3_1=((ps[:,None]*3)*UNd*Ntd+(qs[:,None]*3)*UNd*Ud+k3hs[:,None]*(Ud*Ud*Ud)).sum(0)
        X1=(UNd*2)*ps[:,None]+(Ud*Ud)*qs[:,None]
        X2=(Ntd*Ntd)*ps[:,None]+(UNd*2)*qs[:,None]+(Ud*Ud)*k3hs[:,None]
        K21=X1.T@Ntd+X2.T@Ud; np.fill_diagonal(K21,0.0)
        z=Zs[l].astype(np.float64); a=z.mean(0); zc=z-a; k3s=(zc**3).mean(0)
        k3=k3_1+SCHED[l-1,0]*t*s3; g4=SCHED[l-1,1]
        resid=k3s-k3; noise=6*np.mean(s3*s3)/N_MC; sig=max(np.mean(resid**2)-noise,0)
        k3=k3+(sig/(sig+noise))*resid
        k4v=g4*s4
        if ora_k3 is not None and ora_k3(l): k3=lam[3][l]*s3
        if ora_k4 is not None and ora_k4(l): k4v=lam[4][l]*s4
        muh=mu/np.linalg.norm(mu); Ap=zc@muh
        u=0.75*(C2P[l-1]*s3)+0.25*((zc*zc).T@Ap/N_MC-K21@muh)
        K21=K21+np.outer(u,muh); np.fill_diagonal(K21,0.0)
        t2=t*t; he2=t2-1; he4=(t2-6)*t2+3
        dE1=-k3*t*ph/(6*sz2)+(k3*k3)*he4*ph/(72*s4*sz)+k4v*he2*ph/(24*s3)
        if hi_order:
            H=he_all(-t,8)
            dE1=dE1 + sz*ph*( lam[5][l]/120.0*H[3] + (lam[6][l]/720.0)*H[4]
                            + lam[7][l]/5040.0*H[5] + lam[8][l]/40320.0*H[6] )
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
        m=m_next; Sig=Cn; outs.append(m)
    return np.mean((outs[-1]-Yall[k][-1])**2)

LAM={}
for k in NETS:
    fn=f"rw_{k}_{NB}.npz"
    if not os.path.exists(fn): print("missing cache",fn); sys.exit(1)
    d=np.load(fn); Hm=d['Hm']
    lam={3:Hm[:,3],4:Hm[:,4],5:Hm[:,5],6:Hm[:,6]-10*Hm[:,3]**2,
         7:Hm[:,7]-35*Hm[:,3]*Hm[:,4],8:Hm[:,8]-56*Hm[:,3]*Hm[:,5]-35*Hm[:,4]**2+0*Hm[:,3]}
    LAM[k]=lam
ALL=lambda l: True
LATE=lambda l: l>=10
EARLY=lambda l: l<10
CFG=[("deployed",                 dict()),
     ("exact k3, layers >= 10",   dict(ora_k3=LATE)),
     ("exact k3, layers < 10",    dict(ora_k3=EARLY)),
     ("exact k3, all layers",     dict(ora_k3=ALL)),
     ("exact k4, all layers",     dict(ora_k4=ALL)),
     ("exact k3 and k4",          dict(ora_k3=ALL,ora_k4=ALL)),
     ("exact k3,k4 + k5..k8",     dict(ora_k3=ALL,ora_k4=ALL,hi_order=True))]
print(f"Oracle substitution of MEASURED standardised cumulants ({NB} samples).  Final-layer MSE.\n")
print(f"{'configuration':28s}  " + "  ".join(f"net {k}" for k in NETS) + "   mean      / deployed")
base=None
for name,kw in CFG:
    v=[predict(k,lam=LAM[k],**kw) for k in NETS]
    mv=np.mean(v)
    if base is None: base=mv
    print(f"{name:28s}  " + "  ".join(f"{x:.3e}" for x in v) + f"   {mv:.3e}  {mv/base:7.4f}", flush=True)
