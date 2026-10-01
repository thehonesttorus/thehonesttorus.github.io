import numpy as np, sys
sys.path.insert(0,'/home/user/thehonesttorus.github.io/notes/fresh-slate/designs/heisenberg')
from hd import *
def passes(Ws,N,seed,mean=None,chunk=50000):
    L,n,_=Ws.shape; rng=np.random.default_rng(seed); done=0
    acc=dict(m=np.zeros((L,n)),a=np.zeros((L,n)))
    if mean is not None:
        for k in ['c2','s21','s22','s31']: acc[k]=np.zeros((L,n,n))
        acc['d3']=np.zeros((L,n)); acc['d4']=np.zeros((L,n))
    while done<N:
        b=min(chunk,N-done); a=rng.standard_normal((b,n))
        for l in range(L):
            z=a@Ws[l]; a=np.maximum(z,0); acc['a'][l]+=a.sum(0); acc['m'][l]+=z.sum(0)
            if mean is not None:
                x=z-mean[l]; x2=x*x
                acc['c2'][l]+=x.T@x; acc['s21'][l]+=x2.T@x; acc['s22'][l]+=x2.T@x2; acc['s31'][l]+=(x2*x).T@x
                acc['d3'][l]+=(x2*x).sum(0); acc['d4'][l]+=(x2*x2).sum(0)
        done+=b
    return {k:v/N for k,v in acc.items()}
def run(Ws,st,mode):
    L,n,_=Ws.shape; m=np.zeros(n); C=Ws[0].T@Ws[0]; res=[]
    for l in range(L):
        if mode.endswith('trueC'): m,C=st['m'][l],st['c2'][l]
        G=Gauss(m,C,K=21)
        c2=st['c2'][l]; v=np.diag(c2)
        D=st['d3'][l]; S=st['s21'][l]; K4d=st['d4'][l]-3*v*v
        K22=st['s22'][l]-np.outer(v,v)-2*c2*c2; K31=st['s31'][l]-3*v[:,None]*c2
        if l==0 or mode.startswith('closure'): dEa,dC=np.zeros(n),np.zeros((n,n))
        elif mode.startswith('first'): dEa,dC=inject(G,D,S)
        else: dEa,dC=inject_full2(G,D,S,K4d,K22,K31)
        Ea=G.Ea+dEa; res.append(Ea)
        if l+1<L: Ca=G.cov_a()+dC; m=Ea@Ws[l+1]; C=Ws[l+1].T@Ca@Ws[l+1]
    return np.array(res)
n=int(sys.argv[1]); N=int(float(sys.argv[2])); L=16
for seed in range(2):
    rng=np.random.default_rng(100+seed); Ws=rng.normal(size=(L,n,n))*np.sqrt(2/n)
    p1=passes(Ws,N,seed+50); st=passes(Ws,N,seed+60,mean=p1['m']); st['m']=p1['m']; tr=p1['a']
    for mode in ['closure','first','full2','first_trueC','full2_trueC']:
        e=((run(Ws,st,mode)-tr)**2).mean(1)
        print(seed,"%-12s final %.3e layers"%(mode,e[-1])," ".join("%.1e"%x for x in e[[1,3,7,11,15]]),flush=True)
