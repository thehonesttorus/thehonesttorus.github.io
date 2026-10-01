"""Oracle test: inject the TRUE per-layer diag k3, (2,1) slice, diag k4 of z_l (from MC) into the
reference chain with the same Stein injection as HD. Bounds what better sources/transport can buy."""
import numpy as np, sys, json
sys.path.insert(0,'/home/user/thehonesttorus.github.io/notes/fresh-slate/designs/heisenberg'); sys.path.insert(0,'/root/hdw')
from hd import *
def mc_stats(Ws,N,seed=0,chunk=100000):
    L,n,_=Ws.shape; rng=np.random.default_rng(seed)
    S1=np.zeros((L,n)); S2=np.zeros((L,n,n)); S3=np.zeros((L,n,n)); M3=np.zeros((L,n)); M4=np.zeros((L,n)); A1=np.zeros((L,n)); done=0
    while done<N:
        b=min(chunk,N-done); a=rng.standard_normal((b,n))
        for l in range(L):
            z=a@Ws[l]; S1[l]+=z.sum(0); S2[l]+=z.T@z; z2=z*z; S3[l]+=z2.T@z; M4[l]+=(z2*z2).sum(0)
            a=np.maximum(z,0); A1[l]+=a.sum(0)
        done+=b
    m=S1/N; E2=S2/N; E3=S3/N
    C=E2-m[:,:,None]*m[:,None,:]
    out=[]
    for l in range(L):
        mm=m[l]; e2=E2[l]; e3=E3[l]; d2=np.diag(e2)
        # k3(p,p,q) = E[zp^2 zq] - E[zp^2]m_q - 2 E[zp zq] m_p + 2 m_p^2 m_q
        S=e3-d2[:,None]*mm[None,:]-2*e2*mm[:,None]+2*(mm**2)[:,None]*mm[None,:]
        D=np.diag(S).copy()
        e4=M4[l]/N; e3d=np.diag(e3)
        k4=e4-4*mm*e3d+6*mm**2*d2-3*mm**4-3*(d2-mm**2)**2
        out.append((mm,C[l],D,S,k4))
    return out, A1/N
def run(Ws,st,use_true_mc=False,second=True,inj=True):
    L,n,_=Ws.shape; m=np.zeros(n); C=Ws[0].T@Ws[0]; res=[]
    for l in range(L):
        if use_true_mc: m,C=st[l][0],st[l][1]
        G=Gauss(m,C,K=8)
        _,_,D,S,k4=st[l]
        if l>0 and inj:
            dEa,dC=inject(G,D,S)
            if second:
                d2,d2a=inject2(G,D,k4); dEa=dEa+d2; dC[np.diag_indices(n)]+=d2a-2*G.Ea*d2
        else: dEa,dC=np.zeros(n),np.zeros((n,n))
        Ea=G.Ea+dEa; res.append(Ea)
        if l+1<L:
            Ca=G.cov_a()+dC; m=Ea@Ws[l+1]; C=Ws[l+1].T@Ca@Ws[l+1]
    return np.array(res)
n=int(sys.argv[1]); N=int(float(sys.argv[2])); L=16
for seed in range(2):
    rng=np.random.default_rng(100+seed); Ws=rng.normal(size=(L,n,n))*np.sqrt(2/n)
    st,tr=mc_stats(Ws,N,seed=seed+50)
    for name,kw in [("closure",dict(inj=False)),("oracle k3 1st",dict(second=False)),("oracle k3+k4 2nd",dict(second=True)),
                    ("oracle m,C only",dict(use_true_mc=True,inj=False)),("oracle m,C + k3,k4",dict(use_true_mc=True))]:
        e=((run(Ws,st,**kw)-tr)**2).mean(1)
        print(seed,"%-20s final %.3e  layers"%(name,e[-1])," ".join("%.1e"%x for x in e[[1,3,7,11,15]]),flush=True)
