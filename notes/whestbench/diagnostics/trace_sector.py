# How much of our actual third cumulant lives in the Sigma-metric trace sector R(v,Sigma) that the
# homogeneity identity reduces to a rank-two covariance tangent?
# In whitened coordinates w = Sigma^{-1/2}(h - m):
#   b_i = E[w_i ||w||^2]                       (the metric trace of kappa_3)
#   ||kappa_3||^2 = E[(w . w')^3], w,w' iid    (exact for centered w: third moment = third cumulant)
#   ||R(v,I)||^2 = 3(n+2)||v||^2 with v = b/(n+2), so the trace share is 3||b||^2 / ((n+2)||kappa_3||^2)
import numpy as np, time
import twin2 as T
n=1024; L=16
seeds,Yall=T.load_nets(); W=[w.astype(np.float64) for w in T.regen(seeds[0])]
HS=np.load('../pass2_stats.npz')['HS']; Hm=np.load('../pass2_stats.npz')['Hm']
NS=60_000; CH=6000
rng=np.random.default_rng(99)
print("layer | eff rank | trace-sector share of ||kappa_3||^2 | ||b||/sqrt(n) | corr(b, m_white)")
for LAYER in (2,4,6,8,10,12,14,15):
    Sig=HS[LAYER]; mu=Hm[LAYER]
    evals,V=np.linalg.eigh(Sig)
    keep=evals>1e-6*evals.max(); r=int(keep.sum())
    Vk=V[:,keep]; wk=evals[keep]
    Wt=Vk/np.sqrt(wk)                      # whitening map: w = Wt^T (h - mu), dimension r
    acc_b=np.zeros(r); Ws=[]
    N=0
    for _ in range(NS//CH):
        h=rng.standard_normal((CH,n))
        for k in range(LAYER): h=np.maximum(h@W[k],0.0)
        w=(h-mu)@Wt
        acc_b+=(w*(w*w).sum(1,keepdims=True)).sum(0); N+=CH
        if len(Ws)*CH<12000: Ws.append(w)
    b=acc_b/N
    Wc=np.concatenate(Ws,0)[:8000]
    G=Wc@Wc.T; np.fill_diagonal(G,0.0)
    M=Wc.shape[0]
    k3sq=float((G**3).sum()/(M*(M-1)))
    share=3*float(b@b)/((r+2)*k3sq)
    mw=(np.zeros(r))  # whitened mean is zero by construction; compare b against the whitened mean DIRECTION of mu
    muw=Wt.T@mu
    print(f"{LAYER:5d} |   {r:4d}   |            {share:.4f}                |   {np.linalg.norm(b)/np.sqrt(r):.4f}     |   {abs(np.corrcoef(b,muw)[0,1]):.4f}", flush=True)
