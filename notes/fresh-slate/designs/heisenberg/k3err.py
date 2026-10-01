import numpy as np, sys
sys.path.insert(0,'/root/hdw'); sys.path.insert(0,'/home/user/thehonesttorus.github.io/notes/fresh-slate/designs/heisenberg')
from oracle2 import passes
from hd import hd
n=64; L=16; rng=np.random.default_rng(500); Ws=rng.normal(size=(L,n,n))*np.sqrt(2/n)
p1=passes(Ws,3_000_000,50); st=passes(Ws,3_000_000,60,mean=p1['m'])
tr=[]; hd(Ws,A=None,Kt=3,diagrams="star",trace=tr)
off=~np.eye(n,dtype=bool)
for l,D,S in tr:
    if l in (1,2,4,8,12,15):
        Dt=st['d3'][l]; St=st['s21'][l]
        print(l,"D rel err %.2f corr %.3f | S rel err %.2f corr %.3f"%(np.linalg.norm(D-Dt)/np.linalg.norm(Dt),np.corrcoef(D,Dt)[0,1],np.linalg.norm((S-St)[off])/np.linalg.norm(St[off]),np.corrcoef(S[off],St[off])[0,1]))
