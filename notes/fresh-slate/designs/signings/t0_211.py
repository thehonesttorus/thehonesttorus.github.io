import numpy as np, copula
n=12; N=int(2e6)
rng = np.random.default_rng(1)
W = rng.standard_normal((2, n, n)) * np.sqrt(2 / n)
C1 = W[0].T @ W[0]; s = np.sqrt(np.diag(C1)); R = C1 / np.outer(s, s)
m0 = s/np.sqrt(2*np.pi)
pr = copula.profiles(np.zeros(n), s, np.tile([1.,0,0],(n,1)), rmax=2)
h, F2 = pr; mu=h[:,0]
A2 = F2 - 2*mu[:,None]*h; A2[:,0] += mu**2
F = copula.FACT
# Mehler K211 tensor: centre a + centres b,c (paths), no triangles
K = np.zeros((n,n,n))
R0 = R.copy(); np.fill_diagonal(R0,0)
for p in range(1,4):
    for q in range(1,4):
        K += np.einsum('a,ab,ac,b,c->abc', A2[:,p+q]-2*h[:,p]*h[:,q], R0**p, R0**q, h[:,p], h[:,q])/(F[p]*F[q])
        K += np.einsum('a,ab,bc,b,c->abc', A2[:,p], R0**p, R0**q, h[:,p+q], h[:,q])/(F[p]*F[q])
        K += np.einsum('a,ac,cb,c,b->abc', A2[:,p], R0**p, R0**q, h[:,p+q], h[:,q])/(F[p]*F[q])
E2=np.zeros((n,n)); E211=np.zeros((n,n,n)); done=0
while done<N:
    x=rng.standard_normal((1<<17,n)); a=np.maximum(x@W[0],0)-m0
    E2+=a.T@a; E211+=np.einsum('ia,ib,ic->abc',a**2,a,a); done+=1<<17
E2/=done; E211/=done; v=np.diag(E2)
K211 = E211 - v[:,None,None]*E2[None] - 2*E2[:,:,None]*E2[:,None,:]
idx=[(0,1,2),(3,4,5),(1,0,7),(5,9,2),(0,3,3)]
for t in idx: print(t, K[t], K211[t])
def full3(i,j,k,Fa,Fb,Fc,Dm=8):
    tot=0
    for mab in range(Dm+1):
      for mbc in range(Dm+1):
        for mca in range(Dm+1):
          da,db,dc=mab+mca,mab+mbc,mbc+mca
          if max(da,db,dc)>8 or db==0 or dc==0: continue
          tot+=Fa[da]*Fb[db]*Fc[dc]*R[i,j]**mab*R[j,k]**mbc*R[i,k]**mca/(F[mab]*F[mbc]*F[mca])
    return tot
def cov(i,j): return sum(h[i,d]*h[j,d]*R[i,j]**d/F[d] for d in range(1,9))
for t in idx[:4]:
    i,j,k=t
    full = full3(i,j,k,A2[i],h[j],h[k]) - A2[i,0]*cov(j,k) - 2*cov(i,j)*cov(i,k)
    print(t, "full Mehler", full, "MC", K211[t])
