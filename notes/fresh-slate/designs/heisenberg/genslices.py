import numpy as np
def gen_slices(m,C,W,N,seed,chunk=50000):
    n=len(m); rng=np.random.default_rng(seed); Lc=np.linalg.cholesky(C+1e-10*np.eye(n)); done=0
    S1=np.zeros(n); 
    ys=[]
    acc={k:np.zeros((n,n)) for k in ['c2','s21','s22','s31']}; d3=np.zeros(n); d4=np.zeros(n)
    # two passes: mean then centred moments (same rng stream re-seeded)
    for p in range(2):
        rng=np.random.default_rng(seed); done=0
        while done<N:
            b=min(chunk,N-done); z=m+rng.standard_normal((b,n))@Lc.T; y=np.maximum(z,0)@W; done+=b
            if p==0: S1+=y.sum(0); continue
            x=y-mu; x2=x*x
            acc['c2']+=x.T@x; acc['s21']+=x2.T@x; acc['s22']+=x2.T@x2; acc['s31']+=(x2*x).T@x; d3+=(x2*x).sum(0); d4+=(x2*x2).sum(0)
        if p==0: mu=S1/N
    c2=acc['c2']/N; v=np.diag(c2)
    return dict(D=d3/N, S=acc['s21']/N, K4d=d4/N-3*v*v, K22=acc['s22']/N-np.outer(v,v)-2*c2*c2, K31=acc['s31']/N-3*v[:,None]*c2)
