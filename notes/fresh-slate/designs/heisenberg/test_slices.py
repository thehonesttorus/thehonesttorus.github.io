import numpy as np, sys
sys.path.insert(0,'/root/hdw'); sys.path.insert(0,'/home/user/thehonesttorus.github.io/notes/fresh-slate/designs/heisenberg')
from genslices import gen_slices
from hd import *
rng=np.random.default_rng(3); n=64; L=8
Ws=rng.normal(size=(L,n,n))*np.sqrt(2/n)
m=np.zeros(n); C=Ws[0].T@Ws[0]
for l in range(5):
    G=Gauss(m,C,K=8); m=G.Ea@Ws[l+1]; C=Ws[l+1].T@G.cov_a()@Ws[l+1]
G=Gauss(m,C,K=8); W=Ws[6]
print("offdiag |rho| mean/max", np.abs(G.rho-np.eye(n)).mean(), np.abs(G.rho-np.eye(n)).max())
K31,K22,K4d=gen_slices(m,C,W,4_000_000,1)['K31'],None,None
g=gen_slices(m,C,W,4_000_000,1)
for cyc in [False]:
    a31,a22,a4=gen_k4_slices(G,W,cycle=cyc)
    off=~np.eye(n,dtype=bool)
    for nm,x,y in [("K31",a31[off],g['K31'][off]),("K22",a22[off],g['K22'][off]),("K4d",a4,g['K4d'])]:
        print(cyc,nm,"corr %.3f  rms true %.3e rms diff %.3e"%(np.corrcoef(x,y)[0,1],np.sqrt((y**2).mean()),np.sqrt(((x-y)**2).mean())))
g2=gen_slices(m,C,W,4_000_000,2)
for k in ['K31','K22']:
    x=g2[k][off]; y=g[k][off]; print("MC-MC",k,np.corrcoef(x,y)[0,1], np.sqrt(((x-y)**2).mean())/np.sqrt(2))
