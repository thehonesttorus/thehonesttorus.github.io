import numpy as np, time
rng=np.random.default_rng(5)
for n in [256,512,1024]:
    L=16; N=20000; t=time.time()
    Ws=[rng.normal(0,np.sqrt(2/n),(n,n)) for _ in range(L)]
    Zs=[]
    for b in range(N//2000):
        U=rng.normal(size=(2000,n))
        for W in Ws[:-1]: U=np.maximum(U@W,0)
        Zs.append(U@Ws[-1])                 # final preactivations z_k(x), shape (2000,n)
    Z=np.concatenate(Zs); Zc=Z-Z.mean(0); sd=Zc.std(0)
    sk=(Zc**3).mean(0)/sd**3; ku=(Zc**4).mean(0)/sd**4-3
    # second-moment of skewness over neurons, debiased by sampling variance 6/N
    k3sq=np.mean(sk**2)-6/N
    print(f"n={n}: mean skew {sk.mean():+.4f}  rms skew (debiased) {np.sqrt(max(k3sq,0)):.4f}   1/n={1/n:.4f}  | mean excess kurt {ku.mean():+.4f}  sd kurt {ku.std():.4f} (sampling sd {np.sqrt(24/N):.4f})  [{time.time()-t:.0f}s]")
