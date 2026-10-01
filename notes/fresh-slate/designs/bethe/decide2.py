"""Deciding experiment: v4 at n = 1024 with MC-true node beliefs (v, k3, k4) injected from layer L0 on."""
import sys; sys.path.insert(0,'../../bench'); import bench, bethe, numpy as np, warnings, time
warnings.filterwarnings('ignore')
S=bench.load_set('w1024_d16'); i=int(sys.argv[1]); N=int(float(sys.argv[2]))
W=bench.weights(S,i).astype(np.float64); T=S['means'][i]; L,n,_=W.shape
W32=W.astype(np.float32); rng=np.random.default_rng(21)
acc=np.zeros((L,4,n)); done=0; B=20000; t0=time.time()
while done<N:
    a=rng.standard_normal((B,n),dtype=np.float32)
    for l in range(L):
        z=a@W32[l]; a=np.maximum(z,0); zz=z.astype(np.float64)
        z2=zz*zz; acc[l,0]+=zz.sum(0); acc[l,1]+=z2.sum(0); acc[l,2]+=(z2*zz).sum(0); acc[l,3]+=(z2*z2).sum(0)
    done+=B
print('MC done %.0fs'%(time.time()-t0), flush=True); np.save('/root/mc_nodes_%d.npy'%i, acc/N)
orc=[]
for l in range(L):
    E1,E2,E3,E4=acc[l]/N
    v=E2-E1**2; k3=E3-3*E2*E1+2*E1**3; k4=E4-4*E3*E1-3*E2**2+12*E2*E1**2-6*E1**4
    orc.append(dict(v=v,k3=k3,k4=k4,m=E1))
noise=S['noise'][i]
sub=lambda keys: [{k:d[k] for k in keys} for d in orc]
for name,keys in (('v',('v',)),('k3',('k3',)),('k4',('k4',)),('v+k3',('v','k3')),('k3+k4',('k3','k4')),('v+k4',('v','k4'))):
    e=bethe.estimate_v4(W,node_oracle=sub(keys),oracle_from=6)
    print('v4 + true', name, 'from layer 6: raw %.3e'%(((e[-1]-T[-1])**2).mean()-noise), flush=True)
