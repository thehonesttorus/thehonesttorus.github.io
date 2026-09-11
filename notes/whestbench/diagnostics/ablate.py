# Budget reallocation with an HONEST refit.  Every configuration gets its own kernel scale refitted on
# networks 0-7 (three scalars alpha,gamma,delta multiplying the shipped a/c2p/g4 schedules), then is
# scored on held-out networks 8-15.  Costs are in units of n^3 = 1.07e9 FLOPs; budget 2^41 = 2048 n^3.
#   MC            N*2*L/n                       (2400 -> 75.0)
#   sandwich      3.0 per layer                 (symmetric einsum = 1.5 x 2n^3)
#   source Nt     2.0 per depth level per layer
#   transport     2.0 per extra depth level per layer
#   source K21    4.0 * n_act/n per depth level per layer
# Deployed = 75 + 15*(3+2+3.2) = 198 -> 9.67% of budget, against the harness's measured 9.79%.
import numpy as np, sys, itertools
from scipy.special import ndtr
import twin2 as T
n=1024; L=16; INV=T.INV; pdf=T.pdf
seeds,Yall=T.load_nets()
A_S=np.array([0.00094,0.00783,0.01134,0.01406,0.01658,0.01853,0.02040,0.02217,0.02331,0.02464,
              0.02555,0.02684,0.02781,0.02838,0.02940])
G_S=np.array([0.00715,0.01340,0.01947,0.02302,0.02679,0.03065,0.03478,0.03824,0.04162,0.04615,
              0.04940,0.05096,0.05280,0.05462,0.05873])
C_S=np.array([0.00396,0.07819,0.14538,0.22329,0.30105,0.38174,0.45643,0.55252,0.63162,0.70963,
              0.80605,0.92641,1.00444,1.09326,1.18779])
T_ACT=2.75
WCACHE={}
def getW(k):
    if k not in WCACHE: WCACHE[k]=T.regen(seeds[k])
    return WCACHE[k]

def cost(cfg):
    c=cfg['N']*2.0*L/n
    for l in range(1,L):
        c+=3.0
        d=cfg['depth'](l)
        if d>0:
            c+=2.0*d + 2.0*(d-1)
            if cfg['k21']: c+=3.2*d
    return c

def predict(k, cfg, seed=0):
    W=getW(k); N=cfg['N']; dep=cfg['depth']; k21x=cfg['k21']
    al,ga,de=cfg['al'],cfg['ga'],cfg['de']
    if N>0:
        rng=np.random.default_rng(seed); X=rng.standard_normal((N,n)).astype(np.float32); Zs=[]; h=X
        for l in range(L): z=h@W[l]; Zs.append(z); h=np.maximum(z,0)
    sources={}
    W0=W[0]; G=W0.T@W0; s2=np.diag(G); s=np.sqrt(s2); m=s*INV; so=np.outer(s,s)
    rho=np.clip(G/so,-1,1); C=so*((np.sqrt(1-rho*rho)+rho*(np.pi-np.arccos(rho)))/(2*np.pi))-np.outer(m,m)
    C[np.diag_indices(n)]=s2*(0.5-1/(2*np.pi)); C=C.astype(np.float32); M=(G*0.5).astype(np.float32); M[np.diag_indices(n)]=0
    sources[0]=(M,(INV/s).astype(np.float32),m.astype(np.float32),(s2*s*(2*INV-1.5*INV+2*INV**3)).astype(np.float32),np.full(n,0.5,np.float32)); Sig=C; m=m.astype(np.float32)
    for l in range(1,L):
        Wl=W[l]; mu=(m@Wl).astype(np.float64); Sz=((Wl.T@Sig)@Wl).astype(np.float64); sz2=np.diag(Sz).copy(); sz=np.sqrt(np.maximum(sz2,1e-30)); t=mu/sz
        ph=pdf(t); Ph=ndtr(t); Psi=mu*Ph+sz*ph; m2=(mu*mu+sz2)*Ph+mu*sz*ph; s3=sz2*sz; s4=sz2*sz2
        k3=np.zeros(n,np.float32); K21=np.zeros((n,n),np.float32); d=dep(l)
        if d>0:
            U=Wl; lo=max(l-d,0)
            for srcl in range(l-1,lo-1,-1):
                if srcl not in sources: break
                Ms,ps,qs,k3hs,d1s=sources[srcl]; Nt=Ms.T@U; UN=U*Nt
                k3=k3+((ps[:,None]*np.float32(3))*UN*Nt+(qs[:,None]*np.float32(3))*UN*U+k3hs[:,None]*(U*U*U)).sum(0)
                if k21x:
                    X1=(UN*2)*ps[:,None]+(U*U)*qs[:,None]
                    X2=(Nt*Nt)*ps[:,None]+(UN*2)*qs[:,None]+(U*U)*k3hs[:,None]
                    Kc=X1.T@Nt+X2.T@U; np.fill_diagonal(Kc,0.0); K21+=Kc
                if srcl>lo: U=W[srcl]@(d1s[:,None]*U)
        k3=k3.astype(np.float64)+al*A_S[l-1]*t*s3; g4=de*G_S[l-1]
        muh=mu/np.linalg.norm(mu)
        if N>0:
            z=Zs[l]; a=z.mean(0); zc=z-a; k3s=(zc**3).mean(0)
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
        so=np.outer(sz,sz); rhoc=Sz/so; ds=T.hermite_d(t,ph,Ph,8); rk=rhoc.copy(); Cn=np.zeros((n,n),np.float32); fact=1.0
        for kk in range(1,9):
            fact*=kk; Cn+=rk*np.outer(ds[kk-1],ds[kk-1])/fact
            if kk<8: rk*=rhoc
        Cn*=so; F2=ph/sz
        Cn=Cn+0.5*(K21*np.outer(F2,Ph)+K21.T*np.outer(Ph,F2))
        Cn[np.diag_indices(n)]=(m2+dE2)-m_next*m_next
        Msrc=Sz*Ph[:,None]; np.fill_diagonal(Msrc,0.0); m3=mu*m2+2*sz2*Psi
        sources[l]=(Msrc.astype(np.float32),(ph/sz).astype(np.float32),(2*Psi*(1-Ph)).astype(np.float32),(m3-3*m2*Psi+2*Psi**3).astype(np.float32),Ph.astype(np.float32))
        m=m_next.astype(np.float32); Sig=Cn.astype(np.float32)
    return np.mean((m-Yall[k][-1])**2)

def score(cfg,nets): return float(np.mean([predict(k,cfg) for k in nets]))
FIT=[0,3,6]; HOLD=list(range(8,16))

def refit(cfg):
    best=dict(cfg); bs=score(best,FIT)
    for key,grid in (('al',(0.55,0.75,1.3,1.7)),('ga',(0.7,1.4)),('de',(0.6,1.5))):
        cur=best[key]
        for g in grid:
            c=dict(best); c[key]=cur*g; s=score(c,FIT)
            if s<bs: bs=s; best=c
    return best,bs

CFG={}
def mk(name,depth,k21,N):
    CFG[name]=dict(depth=depth,k21=k21,N=N,al=1.0,ga=1.0,de=1.0)
mk("A deployed",                 lambda l:1, True,  2400)
mk("B no exact K21",             lambda l:1, False, 2400)
mk("C no K21, depth2 from 8",    lambda l:1 if l<8 else 2, False, 2400)
mk("D no K21, depth2 all",       lambda l:2, False, 2400)
mk("E no K21, no MC, depth2 all",lambda l:2, False, 0)
mk("F no K21, MC1200, depth3 all",lambda l:3,False, 1200)
mk("G deployed + depth2 from 12",lambda l:1 if l<12 else 2, True, 2400)
mk("H no K21, depth3 from 6",    lambda l:1 if l<6 else 3, False, 2400)
mk("I deployed, MC 4800",        lambda l:1, True,  4800)
mk("J no K21, no MC, depth3 all",lambda l:3, False, 0)
if __name__=="__main__":
    print(f"{'configuration':32s} {'cost n^3':>9} {'C/B':>7}  {'al':>5} {'ga':>5} {'de':>5}  "
          f"{'fit MSE':>10} {'held MSE':>10} {'adjusted':>10}")
    for name,cfg in CFG.items():
        c=cost(cfg); frac=c/2048.0
        b,bs=refit(cfg)
        hs=score(b,HOLD); adj=hs*max(0.1,frac)
        print(f"{name:32s} {c:9.1f} {frac*100:6.2f}%  {b['al']:5.2f} {b['ga']:5.2f} {b['de']:5.2f}  "
              f"{bs:10.3e} {hs:10.3e} {adj:10.3e}", flush=True)
    