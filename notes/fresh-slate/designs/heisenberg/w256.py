import numpy as np, sys, time, json
sys.path.insert(0,'/root/hdw'); sys.path.insert(0,'/home/user/thehonesttorus.github.io/notes/fresh-slate/designs/heisenberg')
from hd import *; from mc import mc_truth
n=256; L=16
for seed in range(int(sys.argv[1]),int(sys.argv[2])):
    rng=np.random.default_rng(500+seed); Ws=rng.normal(size=(L,n,n))*np.sqrt(2/n)
    tr,noise=mc_truth(Ws,int(2e6),seed=seed)
    cl=closure(Ws); t=time.time(); h=hd(Ws,A=None,Kt=3); dt=time.time()-t
    ec=((cl-tr)**2).mean(1); eh=((h-tr)**2).mean(1)
    print(json.dumps(dict(seed=seed,noise=noise,closure=ec[-1],hd=eh[-1],sec=dt,cl_layers=list(ec),hd_layers=list(eh))),flush=True)
