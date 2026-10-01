import numpy as np, sys, time
sys.path.insert(0,'/home/user/thehonesttorus.github.io/notes/fresh-slate/designs/heisenberg'); sys.path.insert(0,'/root/hdw')
from hd import *; from mc import mc_truth
n=int(sys.argv[1]); L=16; N=int(float(sys.argv[2]))
for seed in range(2):
    rng=np.random.default_rng(100+seed)
    Ws=rng.normal(size=(L,n,n))*np.sqrt(2/n)
    tr,noise=mc_truth(Ws,N,seed=seed)
    cl=closure(Ws)
    res={'closure':cl}
    for A in [0,1,None]:
        t=time.time(); res[f'hd A={A}']=hd(Ws,A=A); 
    for k,v in res.items():
        e=((v-tr)**2).mean(1)
        print(seed,k, "final %.3e"%e[-1], "noise %.1e"%noise, "layers", " ".join("%.1e"%x for x in e[[1,3,7,11,15]]))
