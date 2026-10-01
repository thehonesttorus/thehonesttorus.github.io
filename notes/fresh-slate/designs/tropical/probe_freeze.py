import numpy as np
rng=np.random.default_rng(0)
n=256;L=16;N=40000
W=[rng.standard_normal((n,n))*np.sqrt(2/n) for _ in range(L)]
a=rng.standard_normal((N,n))
for l in range(L):
    z=a@W[l]; a=np.maximum(z,0)
    mu=z.mean(0); s=z.std(0); t=np.abs(mu)/s
    m=a.mean(0); rho=(m**2).sum()/(a**2).mean(0).sum()
    print(l+1, "rho=%.3f"%rho, "hot(|mu/s|<1)=%.2f"%(t<1).mean(), "frozen(|mu/s|>3)=%.2f"%(t>3).mean(), "pon mean %.2f"%(z>0).mean())
