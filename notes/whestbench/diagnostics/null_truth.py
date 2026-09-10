# Same v, same reference: bound the TRUE response by mixed finite differences with common random
# numbers, which is far quieter than importance weighting.  delta_{R(v,Sigma)} M = D_m[v] D_Sigma[Sigma] M.
import numpy as np, time
import twin2 as T
n=1024; L=16
seeds,Yall=T.load_nets(); W=[w.astype(np.float64) for w in T.regen(seeds[0])]
HS=np.load('../pass2_stats.npz')['HS']
def psqrt(S):
    ev,V=np.linalg.eigh(S); ev=np.clip(ev,0,None); return (V*np.sqrt(ev))@V.T
LAYER=8
Sig=HS[LAYER].copy(); Sig=Sig+1e-3*np.mean(np.diag(Sig))*np.eye(n)
rs=np.random.default_rng(3); v=rs.standard_normal(n); v/=np.linalg.norm(v)   # SAME v as the closure test
def suffix(H,l):
    for k in range(l,L): H=np.maximum(H@W[k],0.0)
    return H.mean(0)
NS=200_000; CH=20_000
rng=np.random.default_rng(2718); XI=[rng.standard_normal((CH,n)) for _ in range(NS//CH)]
def Mval(mm,SS):
    A=psqrt(SS); acc=np.zeros(n)
    for xi in XI: acc+=suffix(mm[None,:]+xi@A,LAYER)
    return acc/len(XI)
m0=np.zeros(n)
base=Mval(m0,Sig)
print(f"centered reference, layer {LAYER}, ||M|| = {np.linalg.norm(base):.4e}, {NS} samples, common random numbers")
for s in (0.05,0.10,0.20):
    t=s
    mixed=(Mval(m0+s*v,(1+t)*Sig)-Mval(m0+s*v,(1-t)*Sig)-Mval(m0-s*v,(1+t)*Sig)+Mval(m0-s*v,(1-t)*Sig))/(4*s*t)
    print(f"  s=t={s:.2f}: ||true delta_R M|| = {np.linalg.norm(mixed):.4e}   relative to ||M|| = {np.linalg.norm(mixed)/np.linalg.norm(base):.2e}")
print(f"\n  closure response for the same injection = 8.5086e-02, relative to ||M|| = 6.53e-03")
