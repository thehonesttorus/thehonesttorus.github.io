"""Is the joint kappa_4 that matters local (age 0)?  Replace the true K22/K31 by the slices generated in ONE step
from a Gaussian with the true (m, C) of z_{k-1} (MC)."""
import numpy as np, sys
sys.path.insert(0,'/root/hdw'); sys.path.insert(0,'/home/user/thehonesttorus.github.io/notes/fresh-slate/designs/heisenberg')
from oracle2 import passes
from hd import *
def gen_slices(m,C,W,N,seed,chunk=50000):
    n=len(m); rng=np.random.default_rng(seed); Lc=np.linalg.cholesky(C+1e-10*np.eye(n)); done=0
    S1=np.zeros(n); 
    ys=[]
    acc={k:np.zeros((n,n)) for k in ['c2','s21','s22','s31']}; d3=np.zeros(n); d4=np.zeros(n)
    # two passes: mean then centred moments (same rng stream re-seeded)
    for p in range(2):
        rng=np.random.default_rng(seed); done=0
        while done<N:
            b=min(chunk,N-done); z=m+rng.standard_normal((b,n))@Lc.T; y=np.maximum(z,0)@W; done+=b
            if p==0: S1+=y.sum(0); continue
            x=y-mu; x2=x*x
            acc['c2']+=x.T@x; acc['s21']+=x2.T@x; acc['s22']+=x2.T@x2; acc['s31']+=(x2*x).T@x; d3+=(x2*x).sum(0); d4+=(x2*x2).sum(0)
        if p==0: mu=S1/N
    c2=acc['c2']/N; v=np.diag(c2)
    return dict(D=d3/N, S=acc['s21']/N, K4d=d4/N-3*v*v, K22=acc['s22']/N-np.outer(v,v)-2*c2*c2, K31=acc['s31']/N-3*v[:,None]*c2)
n=int(sys.argv[1]); N=int(float(sys.argv[2])); L=16
for seed in range(2):
    rng=np.random.default_rng(100+seed); Ws=rng.normal(size=(L,n,n))*np.sqrt(2/n)
    p1=passes(Ws,N,seed+50); st=passes(Ws,N,seed+60,mean=p1['m']); st['m']=p1['m']; tr=p1['a']
    gen=[None]+[gen_slices(st['m'][l-1],st['c2'][l-1],Ws[l],N,seed+70+l) for l in range(1,L)]
    for mode in ['full2_true','full2_genJ','full2_genAll']:
        m=np.zeros(n); C=Ws[0].T@Ws[0]; res=[]
        for l in range(L):
            G=Gauss(m,C,K=21); c2=st['c2'][l]; v=np.diag(c2)
            T=dict(D=st['d3'][l],S=st['s21'][l],K4d=st['d4'][l]-3*v*v,K22=st['s22'][l]-np.outer(v,v)-2*c2*c2,K31=st['s31'][l]-3*v[:,None]*c2)
            if l==0: dEa,dC=np.zeros(n),np.zeros((n,n))
            else:
                g=gen[l]
                if mode=='full2_genJ': T['K22'],T['K31']=g['K22'],g['K31']
                if mode=='full2_genAll': T=g
                dEa,dC=inject_full2(G,T['D'],T['S'],T['K4d'],T['K22'],T['K31'])
            Ea=G.Ea+dEa; res.append(Ea)
            if l+1<L: Ca=G.cov_a()+dC; m=Ea@Ws[l+1]; C=Ws[l+1].T@Ca@Ws[l+1]
        e=((np.array(res)-tr)**2).mean(1)
        print(seed,"%-12s final %.3e layers"%(mode,e[-1])," ".join("%.1e"%x for x in e[[1,3,7,11,15]]),flush=True)
