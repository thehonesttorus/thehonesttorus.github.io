import numpy as np
def mc_truth(Ws, N, seed=0, chunk=200000):
    L,n,_=Ws.shape; rng=np.random.default_rng(seed)
    W32=Ws.astype(np.float32); S=np.zeros((L,n)); S2=np.zeros(n); done=0
    while done<N:
        b=min(chunk,N-done); a=rng.standard_normal((b,n),dtype=np.float32)
        for l in range(L):
            a=np.maximum(a@W32[l],0); S[l]+=a.sum(0,dtype=np.float64)
        S2+=(a.astype(np.float64)**2).sum(0); done+=b
    mean=S/N; var=S2/N-mean[-1]**2
    return mean, var.mean()/N
