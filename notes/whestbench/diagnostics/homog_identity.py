# Test the proposed homogeneity identity through a real bias-free ReLU suffix.
#   For f positively homogeneous of degree 1 (every bias-free ReLU suffix coordinate is):
#     D_m[v] D_Sigma[Sigma] M  =  -(1/2) D^2_m[v,m] M  =  D_Sigma[ -(1/2)(v m^T + m v^T) ] M
#   so the covariance-matched trace-cubic response R(v,Sigma) is a rank-<=2 covariance tangent,
#   and vanishes at a centered reference m = 0.
# Central finite differences on Monte Carlo with common random numbers, float64.
import numpy as np, time
import twin2 as T
n=1024; L=16
seeds,Yall=T.load_nets(); W=[w.astype(np.float64) for w in T.regen(seeds[0])]; Y=Yall[0].astype(np.float64)
HS=np.load('../pass2_stats.npz')['HS']
def suffix(H,l):
    for k in range(l,L): H=np.maximum(H@W[k],0.0)
    return H.mean(0)
def psqrt(S):
    w,V=np.linalg.eigh(S); w=np.clip(w,0.0,None); return (V*np.sqrt(w))@V.T
LAYER=8
m=Y[LAYER-1].copy(); Sig=HS[LAYER].copy()
lam=1e-3*np.mean(np.diag(Sig)); Sig=Sig+lam*np.eye(n)
ev=np.linalg.eigvalsh(Sig)
print(f"layer {LAYER}: ridged Sigma eig min {ev[0]:.3e} max {ev[-1]:.3e}; ||m|| {np.linalg.norm(m):.3f}")
NS=200_000; CH=20_000
rng=np.random.default_rng(4242)
XI=[rng.standard_normal((CH,n)) for _ in range(NS//CH)]
def Mval(mm,SS):
    A=psqrt(SS); acc=np.zeros(n)
    for xi in XI: acc+=suffix(mm[None,:]+xi@A, LAYER)
    return acc/len(XI)
rs=np.random.default_rng(7); v=rs.standard_normal(n); v/=np.linalg.norm(v)
s=t=0.05
t0=time.time()
mixed=(Mval(m+s*v,(1+t)*Sig)-Mval(m+s*v,(1-t)*Sig)-Mval(m-s*v,(1+t)*Sig)+Mval(m-s*v,(1-t)*Sig))/(4*s*t)
dS=-0.5*(np.outer(v,m)+np.outer(m,v))
covtan=(Mval(m,Sig+t*dS)-Mval(m,Sig-t*dS))/(2*t)
d2=(Mval(m+s*v+s*m,Sig)-Mval(m+s*v-s*m,Sig)-Mval(m-s*v+s*m,Sig)+Mval(m-s*v-s*m,Sig))/(4*s*s)
def cmp(a,b,na,nb):
    print(f"  {na:26s} vs {nb:26s}: rel diff {np.linalg.norm(a-b)/(0.5*(np.linalg.norm(a)+np.linalg.norm(b))):.4f}"
          f"   corr {np.corrcoef(a,b)[0,1]:.5f}   norms {np.linalg.norm(a):.4e} / {np.linalg.norm(b):.4e}")
print(f"[{time.time()-t0:.0f}s] suffix h_{LAYER} -> h_16, {NS} samples, common random numbers, s=t={s}")
cmp(mixed,covtan,"D_m[v]D_Sig[Sig]M","D_Sig[-(vm^T+mv^T)/2]M")
cmp(mixed,-0.5*d2,"D_m[v]D_Sig[Sig]M","-(1/2)D^2_m[v,m]M")
m0=np.zeros(n)
mixed0=(Mval(m0+s*v,(1+t)*Sig)-Mval(m0+s*v,(1-t)*Sig)-Mval(m0-s*v,(1+t)*Sig)+Mval(m0-s*v,(1-t)*Sig))/(4*s*t)
ref=Mval(m0,Sig)
print(f"  centered reference m=0: ||mixed|| {np.linalg.norm(mixed0):.4e} vs ||M(0,Sigma)|| {np.linalg.norm(ref):.4e}"
      f"  ratio {np.linalg.norm(mixed0)/np.linalg.norm(ref):.2e}   (identity predicts 0)")
print(f"  for scale, at m != 0: ||mixed||/||M(m,Sigma)|| = {np.linalg.norm(mixed)/np.linalg.norm(Mval(m,Sig)):.2e}")
