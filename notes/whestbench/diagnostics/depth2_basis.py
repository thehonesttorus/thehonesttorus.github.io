# The depth-2 source increment: what IS it, and does it lie in a low-dimensional analytic basis?
# Runs the DEPLOYED trajectory (depth 1 + shipped kernels) and at each layer additionally forms the
# source term that a depth-2 channel would add.  Then asks whether that increment is captured by
# sigma^3 times a low-order polynomial in t -- and, since the shipped kernel already covers t*sigma^3,
# what is left after projecting that direction out.
import numpy as np, sys
from scipy.special import ndtr
import twin2 as T
n=1024; L=16; INV=T.INV; pdf=T.pdf
seeds,Yall=T.load_nets()
SCHED=np.array([[0.00094,0.00715],[0.00783,0.01340],[0.01134,0.01947],[0.01406,0.02302],[0.01658,0.02679],
 [0.01853,0.03065],[0.02040,0.03478],[0.02217,0.03824],[0.02331,0.04162],[0.02464,0.04615],
 [0.02555,0.04940],[0.02684,0.05096],[0.02781,0.05280],[0.02838,0.05462],[0.02940,0.05873]])
C2P=np.array([0.00396,0.07819,0.14538,0.22329,0.30105,0.38174,0.45643,0.55252,0.63162,0.70963,
              0.80605,0.92641,1.00444,1.18779-0.09,1.18779])
C2P=np.array([0.00396,0.07819,0.14538,0.22329,0.30105,0.38174,0.45643,0.55252,0.63162,0.70963,
              0.80605,0.92641,1.00444,1.09326,1.18779])

def trace(k, N_MC=2400, seed=0):
    """Deployed trajectory; returns per-layer (t, sz, dk3_from_l2, du_from_l2, k3_dep, u_dep)."""
    W=T.regen(seeds[k]); rng=np.random.default_rng(seed)
    X=rng.standard_normal((N_MC,n)).astype(np.float32); Zs=[]; h=X
    for l in range(L): z=h@W[l]; Zs.append(z); h=np.maximum(z,0)
    rec={}; sources={}
    W0=W[0]; G=W0.T@W0; s2=np.diag(G); s=np.sqrt(s2); m=s*INV; so=np.outer(s,s)
    rho=np.clip(G/so,-1,1); C=so*((np.sqrt(1-rho*rho)+rho*(np.pi-np.arccos(rho)))/(2*np.pi))-np.outer(m,m)
    C[np.diag_indices(n)]=s2*(0.5-1/(2*np.pi)); M=G*0.5; M[np.diag_indices(n)]=0
    sources[0]=(M,INV/s,m,s2*s*(2*INV-1.5*INV+2*INV**3),np.full(n,0.5)); Sig=C
    for l in range(1,L):
        Wl=W[l]; mu=m@Wl; Sz=(Wl.T@Sig)@Wl; sz2=np.diag(Sz).copy(); sz=np.sqrt(np.maximum(sz2,1e-30)); t=mu/sz
        ph=pdf(t); Ph=ndtr(t); Psi=mu*Ph+sz*ph; m2=(mu*mu+sz2)*Ph+mu*sz*ph; s3=sz2*sz; s4=sz2*sz2
        # --- depth 1 (deployed) ---
        U=Wl.astype(np.float64); Ms,ps,qs,k3hs,d1s=sources[l-1]
        Nt=Ms.T@U; UN=U*Nt
        k3_1=((ps[:,None]*3)*UN*Nt+(qs[:,None]*3)*UN*U+k3hs[:,None]*(U*U*U)).sum(0)
        X1=(UN*2)*ps[:,None]+(U*U)*qs[:,None]; X2=(Nt*Nt)*ps[:,None]+(UN*2)*qs[:,None]+(U*U)*k3hs[:,None]
        K21_1=X1.T@Nt+X2.T@U; np.fill_diagonal(K21_1,0.0)
        # --- depth 2 increment (not used in the trajectory) ---
        dk3=np.zeros(n); dK21=None
        if l>=2:
            U2=W[l-1]@(d1s[:,None]*U); Ms2,ps2,qs2,k3hs2,d1s2=sources[l-2]
            Nt2=Ms2.T@U2; UN2=U2*Nt2
            dk3=((ps2[:,None]*3)*UN2*Nt2+(qs2[:,None]*3)*UN2*U2+k3hs2[:,None]*(U2*U2*U2)).sum(0)
            Y1=(UN2*2)*ps2[:,None]+(U2*U2)*qs2[:,None]
            Y2=(Nt2*Nt2)*ps2[:,None]+(UN2*2)*qs2[:,None]+(U2*U2)*k3hs2[:,None]
            dK21=Y1.T@Nt2+Y2.T@U2; np.fill_diagonal(dK21,0.0)
        # --- continue on the deployed trajectory ---
        k3=k3_1+SCHED[l-1,0]*t*s3; g4=SCHED[l-1,1]
        z=Zs[l].astype(np.float64); a=z.mean(0); zc=z-a; k3s=(zc**3).mean(0)
        resid=k3s-k3; noise=6*np.mean(s3*s3)/N_MC; sig=max(np.mean(resid**2)-noise,0)
        k3=k3+(sig/(sig+noise))*resid
        muh=mu/np.linalg.norm(mu); Ap=zc@muh
        u_mc=(zc*zc).T@Ap/N_MC-K21_1@muh
        u=0.75*(C2P[l-1]*s3)+0.25*u_mc
        rec[l]=dict(t=t,sz=sz,s3=s3,dk3=dk3,k3=k3.copy(),k3_1=k3_1.copy(),
                    du=(dK21@muh if dK21 is not None else np.zeros(n)),u=u.copy(),
                    k3_true=k3s.copy(),muh=muh)
        K21=K21_1+np.outer(u,muh); np.fill_diagonal(K21,0.0)
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
        Msrc=Sz*Ph[:,None]; np.fill_diagonal(Msrc,0.0)
        m3=mu*m2+2*sz2*Psi
        sources[l]=(Msrc,ph/sz,2*Psi*(1-Ph),m3-3*m2*Psi+2*Psi**3,Ph)
        m=m_next; Sig=Cn
    return rec

def r2(y,X):
    c,_,_,_=np.linalg.lstsq(X,y,rcond=None); r=y-X@c
    return 1-np.sum(r*r)/np.sum(y*y), c

NETS=list(range(8))
recs=[trace(k) for k in NETS]
print("DEPTH-2 SOURCE INCREMENT: basis analysis, networks 0-7\n")
print("A. relative size of the depth-2 increment against the depth-1 term it sits on")
print(f"{'layer':>5} {'||dk3||/||k3_1||':>17} {'||du||/||u||':>14}")
for l in range(2,L):
    a=np.mean([np.linalg.norm(r[l]['dk3'])/np.linalg.norm(r[l]['k3_1']) for r in recs])
    b=np.mean([np.linalg.norm(r[l]['du'])/np.linalg.norm(r[l]['u']) for r in recs])
    print(f"{l:5d} {a:17.4f} {b:14.4f}")

print("\nB. R^2 of the depth-2 k3 increment against analytic bases (mean over nets 0-7)")
print(f"{'layer':>5} {'t*s3':>8} {'s3':>8} {'t,1':>8} {'t,1,t2':>9} {'t,1,t2,t3':>11} {'+phi':>8}")
for l in range(2,L):
    rows=[]
    for r in recs:
        t=r[l]['t']; s3=r[l]['s3']; y=r[l]['dk3']; ph=pdf(t)
        b_t=(t*s3)[:,None]; b_1=s3[:,None]; b_t2=(t*t*s3)[:,None]; b_t3=(t**3*s3)[:,None]; b_p=(ph*s3)[:,None]
        rows.append([r2(y,b_t)[0], r2(y,b_1)[0], r2(y,np.hstack([b_t,b_1]))[0],
                     r2(y,np.hstack([b_t,b_1,b_t2]))[0], r2(y,np.hstack([b_t,b_1,b_t2,b_t3]))[0],
                     r2(y,np.hstack([b_t,b_1,b_t2,b_t3,b_p]))[0]])
    m=np.mean(rows,0)
    print(f"{l:5d} " + " ".join(f"{v:8.4f}" for v in m[:2]) + f" {m[2]:8.4f} {m[3]:9.4f} {m[4]:11.4f} {m[5]:8.4f}")
