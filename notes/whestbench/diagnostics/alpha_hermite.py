# Is alpha > 1 on the second Price correction a self-energy, or is it approximating the exact
# Hermite tail (k >= 3) that a two-term truncation drops?
#   exact ReLU covariance:  C_ij = sigma_i sigma_j * sum_k rho^k d_k(t_i) d_k(t_j) / k!
#   k=1 is the "first order" term, k=2 is the second Price correction, k>=3 is the dropped tail.
# alpha_implied = (sum_{k>=2}) / (k=2 term):  the alpha a two-term closure would need to match exact.
import numpy as np
from scipy.special import ndtr
import twin2 as T
n=1024; L=16; pdf=T.pdf
seeds,Yall=T.load_nets(); HS=np.load('../pass2_stats.npz')['HS']
W=T.regen(seeds[0]); Y=Yall[0].astype(np.float64)
print("Implied alpha on the second Price term, from the exact Hermite series (network 0).")
print("weighted by |k=2 term| so it reflects the pairs that matter; |rho| stats for context.\n")
print("layer | rms|rho| | max|rho| | implied alpha (weighted) | implied alpha at rho=rms | share of tail in k=3 vs k>=4")
for l in range(1,L):
    mu=Y[l-1]@W[l].astype(np.float64); Sz=W[l].T.astype(np.float64)@HS[l]@W[l].astype(np.float64)
    sg=np.sqrt(np.diag(Sz)); t=mu/sg; rho=Sz/np.outer(sg,sg); np.fill_diagonal(rho,0.0)
    ph=pdf(t); Ph=ndtr(t); ds=T.hermite_d(t,ph,Ph,10)
    terms=[]
    rk=rho.copy(); fact=1.0
    for k in range(1,11):
        fact*=k; terms.append(rk*np.outer(ds[k-1],ds[k-1])/fact)
        rk=rk*rho
    t2=terms[1]; tail=sum(terms[1:])           # k=2 alone, and k>=2 total
    wsum=np.abs(t2).sum()
    alpha_w=np.sum(np.abs(t2)*np.where(np.abs(t2)>0,tail/np.where(t2==0,1,t2),0))/max(wsum,1e-30)
    r0=np.sqrt(np.mean(rho**2))
    # implied alpha at a representative correlation, using median t
    tm=np.median(t); phm=pdf(tm); Phm=ndtr(tm); dm=T.hermite_d(np.array([tm]),np.array([phm]),np.array([Phm]),10)
    num=sum((r0**k)*dm[k-1][0]**2/np.math.factorial(k) if hasattr(np,'math') else (r0**k)*dm[k-1][0]**2/float(__import__('math').factorial(k)) for k in range(2,11))
    den=(r0**2)*dm[1][0]**2/2.0
    k3s=np.abs(terms[2]).sum(); k4s=sum(np.abs(x).sum() for x in terms[3:])
    print(f"{l+1:5d} |  {r0:.4f}  |  {np.abs(rho).max():.3f}   |        {alpha_w:.3f}          |        {num/den:.3f}          | {k3s/(k3s+k4s):.2f} / {k4s/(k3s+k4s):.2f}")
