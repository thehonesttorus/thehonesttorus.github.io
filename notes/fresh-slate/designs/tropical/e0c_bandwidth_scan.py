import numpy as np
rng=np.random.default_rng(3); n=16; L=2
W=[rng.standard_normal((n,n))*np.sqrt(2/n) for _ in range(L)]
hs=np.array([0.005,0.01,0.02,0.04,0.08])
lhs=np.zeros(n); rhs=np.zeros((len(hs),n)); N=0
X=rng.standard_normal((100000,n)); s=[X@W[0]]; s.append((np.maximum(s[0],0)@W[1])); s=[v.std(0) for v in s]
for b in range(100):
    X=rng.standard_normal((20000,n)); z1=X@W[0]; g1=(z1>0); a1=z1*g1; z2=a1@W[1]; g2=z2>0
    lhs+=(z2*g2).sum(0); N+=20000
    G2=(g1[:,:,None]*W[1][None])  # d z2/d a1-ish ; grad_x z2 = W0 @ diag(g1) @ W1
    gr2=np.einsum('in,bnj->bij',W[0],G2); n2=(gr2**2).sum(1)
    n1=(W[0]**2).sum(0)
    for hi,h in enumerate(hs):
        k1=np.exp(-.5*(z1/(h*s[0]))**2)/(np.sqrt(2*np.pi)*h*s[0])*n1
        k2=np.exp(-.5*(z2/(h*s[1]))**2)/(np.sqrt(2*np.pi)*h*s[1])*n2
        rhs[hi]+= k2.sum(0) + ((k1@W[1])*g2).sum(0)
lhs/=N; rhs/=N
for hi,h in enumerate(hs): print(h, "mean(RHS-LHS)/mean(LHS) = %.4f  max abs %.4f"%((rhs[hi]-lhs).mean()/lhs.mean(), np.abs(rhs[hi]-lhs).max()))
