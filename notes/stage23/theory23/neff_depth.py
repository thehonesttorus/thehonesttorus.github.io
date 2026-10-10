import numpy as np
rng=np.random.default_rng(7)
n=1024; L=16; N=4000
X=rng.normal(size=(N,n)); U=X
print(" l   PR(C_l full)  PR(C_perp)  n/PR_perp   prediction n/PR ~ l+1 (linear free-MP step)")
for l in range(1,L):
    W=rng.normal(0,np.sqrt(2/n),(n,n)); U=np.maximum(U@W,0)
    mu=U.mean(0); muh=mu/np.linalg.norm(mu)
    Uc=U-mu; Up=Uc-np.outer(Uc@muh,muh)
    # PR via Gram trick (N>n here so use covariance directly)
    C=Up.T@Up/N; tr=np.trace(C); fr=np.sum(C*C)
    Cf=Uc.T@Uc/N; prf=np.trace(Cf)**2/np.sum(Cf*Cf)
    # finite-sample bias of ||C||_F^2: subtract tr^2/N approx
    frc=fr-tr**2/N
    print(f"{l:2d}   {prf:8.1f}    {tr**2/frc:8.1f}    {n/(tr**2/frc):6.2f}       {l+1}")
