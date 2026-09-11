# The depth-2 increment is ~90% the size of the depth-1 term and only 25% captured by the (t,sigma)
# basis.  So measure the WHOLE series: at a deep layer, the source contribution from every earlier
# layer, its norm, the running partial sum, the coherence, and how well each partial sum is captured
# by the analytic basis the shipped kernel uses.
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
TARGET=[int(x) for x in sys.argv[1].split(',')] if len(sys.argv)>1 else [8,12,15]

def run(k, N_MC=2400, seed=0):
    W=T.regen(seeds[k]); rng=np.random.default_rng(seed)
    X=rng.standard_normal((N_MC,n)).astype(np.float32); Zs=[]; h=X
    for l in range(L): z=h@W[l]; Zs.append(z); h=np.maximum(z,0)
    out={}; sources={}
    W0=W[0]; G=W0.T@W0; s2=np.diag(G); s=np.sqrt(s2); m=s*INV; so=np.outer(s,s)
    rho=np.clip(G/so,-1,1); C=so*((np.sqrt(1-rho*rho)+rho*(np.pi-np.arccos(rho)))/(2*np.pi))-np.outer(m,m)
    C[np.diag_indices(n)]=s2*(0.5-1/(2*np.pi)); M=G*0.5; M[np.diag_indices(n)]=0
    sources[0]=(M,INV/s,m,s2*s*(2*INV-1.5*INV+2*INV**3),np.full(n,0.5)); Sig=C
    for l in range(1,L):
        Wl=W[l]; mu=m@Wl; Sz=(Wl.T@Sig)@Wl; sz2=np.diag(Sz).copy(); sz=np.sqrt(np.maximum(sz2,1e-30)); t=mu/sz
        ph=pdf(t); Ph=ndtr(t); Psi=mu*Ph+sz*ph; m2=(mu*mu+sz2)*Ph+mu*sz*ph; s3=sz2*sz; s4=sz2*sz2
        U=Wl.astype(np.float64); terms=[]
        depth = l if l in TARGET else 1
        for srcl in range(l-1,l-1-depth,-1):
            Ms,ps,qs,k3hs,d1s=sources[srcl]; Nt=Ms.T@U; UN=U*Nt
            terms.append(((ps[:,None]*3)*UN*Nt+(qs[:,None]*3)*UN*U+k3hs[:,None]*(U*U*U)).sum(0))
            if srcl>l-depth: U=W[srcl]@(d1s[:,None]*U)
        k3_1=terms[0]
        # deployed trajectory uses depth 1 + kernel
        X1=(Wl.astype(np.float64)*(sources[l-1][0].T@Wl.astype(np.float64)))*2*sources[l-1][1][:,None] \
           +(Wl.astype(np.float64)**2)*sources[l-1][2][:,None]
        Ud=Wl.astype(np.float64); Ms,ps,qs,k3hs,d1s=sources[l-1]; Ntd=Ms.T@Ud; UNd=Ud*Ntd
        X1=(UNd*2)*ps[:,None]+(Ud*Ud)*qs[:,None]
        X2=(Ntd*Ntd)*ps[:,None]+(UNd*2)*qs[:,None]+(Ud*Ud)*k3hs[:,None]
        K21=X1.T@Ntd+X2.T@Ud; np.fill_diagonal(K21,0.0)
        z=Zs[l].astype(np.float64); a=z.mean(0); zc=z-a; k3s=(zc**3).mean(0)
        if l in TARGET: out[l]=dict(terms=[x.copy() for x in terms],t=t.copy(),s3=s3.copy(),truth=k3s.copy())
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

def r2(y,X):
    c,_,_,_=np.linalg.lstsq(X,y,rcond=None); r=y-X@c; return 1-np.sum(r*r)/np.sum(y*y), c

NETS=[0,1,2,3]
R=[run(k) for k in NETS]
for l in TARGET:
    print(f"\n=== layer {l}: the source series, averaged over networks {NETS} ===")
    nt=len(R[0][l]['terms'])
    print(f"{'j (src=l-j)':>11} {'||term_j||':>12} {'/||term_1||':>12} {'||sum_1..j||':>13} {'coh':>7} "
          f"{'R2(term_j|t*s3)':>16} {'R2(sum|t*s3)':>13} {'R2(sum|t,1,t2)':>15}")
    for j in range(1,nt+1):
        tn=[];rel=[];sn=[];coh=[];r2t=[];r2s=[];r2s3=[]
        for rr in R:
            T_=rr[l]['terms']; t=rr[l]['t']; s3=rr[l]['s3']
            S=np.sum(T_[:j],0); n1=np.linalg.norm(T_[0])
            tn.append(np.linalg.norm(T_[j-1])); rel.append(np.linalg.norm(T_[j-1])/n1)
            sn.append(np.linalg.norm(S)/n1)
            coh.append(np.linalg.norm(S)/sum(np.linalg.norm(x) for x in T_[:j]))
            b=np.column_stack([t*s3,s3,t*t*s3])
            r2t.append(r2(T_[j-1],b[:,:1])[0]); r2s.append(r2(S,b[:,:1])[0]); r2s3.append(r2(S,b)[0])
        print(f"{j:11d} {np.mean(tn):12.4e} {np.mean(rel):12.4f} {np.mean(sn):13.4f} {np.mean(coh):7.3f} "
              f"{np.mean(r2t):16.4f} {np.mean(r2s):13.4f} {np.mean(r2s3):15.4f}")
    # how well does the FULL sum match the measured truth, and where does the kernel sit?
    for rr in R[:1]:
        T_=rr[l]['terms']; t=rr[l]['t']; s3=rr[l]['s3']; tr=rr[l]['truth']
        for j in [1,2,4,8,nt]:
            if j>nt: continue
            S=np.sum(T_[:j],0)
            print(f"    depth {j:2d}: ||S-truth||/||truth|| = {np.linalg.norm(S-tr)/np.linalg.norm(tr):.4f}"
                  f"   corr(S,truth) = {np.corrcoef(S,tr)[0,1]:.4f}")
