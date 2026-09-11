# The source series is COHERENT (0.65) and its SUM is 90% rank-one in t*sigma^3.  Two models that
# costs almost nothing follow:
#   (R) renormalisation      k3 ~ gamma_l * (depth-1 term)          [1 multiply/layer]
#   (K) kernel               k3 ~ a_l * t*sigma^3                   [1 multiply/layer, shipped]
#   (RK) both                k3 ~ gamma_l*(depth-1) + a_l*t*sigma^3 [2 multiplies/layer]
# Regress against (i) the exact full-depth source sum and (ii) a high-accuracy MC third cumulant.
import numpy as np, sys
from scipy.special import ndtr
import twin2 as T
n=1024; L=16; INV=T.INV; pdf=T.pdf
seeds,Yall=T.load_nets()
SCHED=np.array([[0.00094,0.00715],[0.00783,0.01340],[0.01134,0.01947],[0.01406,0.02302],[0.01658,0.02679],
 [0.01853,0.03065],[0.02040,0.03478],[0.02217,0.03824],[0.02331,0.04162],[0.02464,0.04615],
 [0.02555,0.04940],[0.02684,0.05096],[0.02781,0.05280],[0.02838,0.05462],[0.02940,0.05873]])
C2P=np.array([0.00396,0.07819,0.14538,0.22329,0.30105,0.38174,0.45643,0.55252,0.63162,0.70963,
              0.80605,0.92641,1.00444,1.09326,1.18779])
NBIG=int(sys.argv[1]) if len(sys.argv)>1 else 60000

def run(k, N_MC=2400, seed=0, NB=NBIG):
    W=T.regen(seeds[k]); rng=np.random.default_rng(seed)
    X=rng.standard_normal((N_MC,n)).astype(np.float32); Zs=[]; h=X
    for l in range(L): z=h@W[l]; Zs.append(z); h=np.maximum(z,0)
    # high-accuracy truth for k3, chunked
    rb=np.random.default_rng(12345+k); S1=np.zeros((L,n)); S2=np.zeros((L,n)); S3=np.zeros((L,n)); ch=6000
    done=0
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
        U=Wl.astype(np.float64); full=np.zeros(n); k3_1=None
        for srcl in range(l-1,-1,-1):
            Ms,ps,qs,k3hs,d1s=sources[srcl]; Nt=Ms.T@U; UN=U*Nt
            tm=((ps[:,None]*3)*UN*Nt+(qs[:,None]*3)*UN*U+k3hs[:,None]*(U*U*U)).sum(0)
            full+=tm
            if k3_1 is None: k3_1=tm.copy()
            if srcl>0: U=W[srcl]@(d1s[:,None]*U)
        out[l]=dict(t=t.copy(),s3=s3.copy(),k3_1=k3_1,full=full,truth=K3T[l].copy())
        # deployed trajectory
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

def fit(y,X):
    c,_,_,_=np.linalg.lstsq(X,y,rcond=None); r=y-X@c; return 1-np.sum(r*r)/np.sum(y*y), c

NETS=[0,1,2,3]
R=[run(k) for k in NETS]
print(f"High-accuracy truth from {NBIG} samples.  Networks {NETS}.\n")
print("A. TARGET = exact full-depth source sum.  R^2 of each cheap model.")
print(f"{'layer':>5} {'K: a*t*s3':>10} {'R: g*k3_1':>10} {'RK: both':>10} {'gamma':>8} {'a_fit':>10} {'a_ship':>10}")
for l in range(1,L):
    r_k=[];r_r=[];r_rk=[];gam=[];af=[]
    for rr in R:
        d=rr[l]; b_k=(d['t']*d['s3'])[:,None]; b_r=d['k3_1'][:,None]; y=d['full']
        r_k.append(fit(y,b_k)[0]); v=fit(y,b_r); r_r.append(v[0]); gam.append(v[1][0])
        v2=fit(y,np.hstack([b_r,b_k])); r_rk.append(v2[0]); af.append(v2[1][1])
    print(f"{l:5d} {np.mean(r_k):10.4f} {np.mean(r_r):10.4f} {np.mean(r_rk):10.4f} {np.mean(gam):8.3f} "
          f"{np.mean(af):10.5f} {SCHED[l-1,0]:10.5f}")

print("\nB. TARGET = measured third cumulant (the thing that actually matters).")
print(f"{'layer':>5} {'depth1 only':>12} {'K only':>9} {'d1+K (ship)':>12} {'R: g*d1':>9} "
      f"{'RK':>8} {'full sum':>10} {'full+K':>9} {'g_fit':>7}")
for l in range(1,L):
    cols=[[] for _ in range(8)]
    for rr in R:
        d=rr[l]; y=d['truth']; b_k=(d['t']*d['s3'])[:,None]; b_1=d['k3_1'][:,None]; b_f=d['full'][:,None]
        # 'depth1 only' and 'd1+K(ship)' use FIXED coefficients (no fitting) -> report 1 - ||y-pred||^2/||y||^2
        def ev(p): r=y-p; return 1-np.sum(r*r)/np.sum(y*y)
        cols[0].append(ev(d['k3_1']))
        cols[1].append(fit(y,b_k)[0])
        cols[2].append(ev(d['k3_1']+SCHED[l-1,0]*d['t']*d['s3']))
        v=fit(y,b_1); cols[3].append(v[0]); cols[7].append(v[1][0])
        cols[4].append(fit(y,np.hstack([b_1,b_k]))[0])
        cols[5].append(ev(d['full']))
        cols[6].append(fit(y,np.hstack([b_f,b_k]))[0])
    m=[np.mean(c) for c in cols]
    print(f"{l:5d} {m[0]:12.4f} {m[1]:9.4f} {m[2]:12.4f} {m[3]:9.4f} {m[4]:8.4f} {m[5]:10.4f} {m[6]:9.4f} {m[7]:7.3f}")
