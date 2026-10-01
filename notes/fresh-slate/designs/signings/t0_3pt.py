import numpy as np, copula, itertools
F=copula.FACT
rho = {('a','b'):0.35, ('b','c'):-0.25, ('a','c'):0.3}
Rm = np.array([[1,.35,.3],[.35,1,-.25],[.3,-.25,1]])
s = np.array([1.2,0.9,1.1])
h, F2 = copula.profiles(np.zeros(3), s, np.tile([1.,0,0],(3,1)), rmax=2, dmax=8)
mu=h[:,0]; A2=F2-2*mu[:,None]*h; A2[:,0]+=mu**2
# Mehler: E[A2(g_a) B(g_b) C(g_c)] with B,C centred
tot=0; Dm=8
for mab in range(0,Dm+1):
  for mbc in range(0,Dm+1):
    for mca in range(0,Dm+1):
      da,db,dc = mab+mca, mab+mbc, mbc+mca
      if max(da,db,dc)>8 or db==0 or dc==0: continue
      tot += A2[0,da]*h[1,db]*h[2,dc]*Rm[0,1]**mab*Rm[1,2]**mbc*Rm[0,2]**mca/(F[mab]*F[mbc]*F[mca])
rng=np.random.default_rng(0); acc=0; N=0; acc2=np.zeros(3)
L=np.linalg.cholesky(Rm)
for _ in range(40):
    g=rng.standard_normal((1<<20,3))@L.T; a=np.maximum(g*s,0)-mu
    acc+=np.sum(a[:,0]**2*a[:,1]*a[:,2]); N+=len(g)
print("Mehler", tot, "MC", acc/N)
