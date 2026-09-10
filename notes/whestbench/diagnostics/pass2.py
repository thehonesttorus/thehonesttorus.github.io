# Accurate per-layer statistics of network 0: means/covariances of every h_l and raw moments (1..4) of every
# pre-activation z_l, from N2 samples. Saved for closure diagnostics (no budget accounting here).
import numpy as np, time
t0=time.time()
W=np.load('data/W0.npy'); n=1024; L=16
rng=np.random.default_rng(2024); ch=8192; N2=1_500_000
Hm=np.zeros((L+1,n)); HS=np.zeros((L+1,n,n)); Sz=np.zeros((L,4,n))
for s in range(0,N2,ch):
    h=rng.standard_normal((ch,n)).astype(np.float32)
    for l in range(L+1):
        hd=h.astype(np.float64); Hm[l]+=hd.sum(0); HS[l]+=hd.T@hd
        if l<L:
            z=h@W[l]; zd=z.astype(np.float64); z2=zd*zd
            Sz[l,0]+=zd.sum(0); Sz[l,1]+=z2.sum(0); Sz[l,2]+=(z2*zd).sum(0); Sz[l,3]+=(z2*z2).sum(0)
            h=np.maximum(z,0)
    if (s//ch)%20==0: print(f"  {s+ch}/{N2} {time.time()-t0:.0f}s", flush=True)
Hm/=N2; HS=HS/N2-np.einsum('li,lj->lij',Hm,Hm); Sz/=N2
np.savez('pass2_stats.npz',Hm=Hm,HS=HS,Sz=Sz,N2=N2)
print(f"done {time.time()-t0:.0f}s", flush=True)
