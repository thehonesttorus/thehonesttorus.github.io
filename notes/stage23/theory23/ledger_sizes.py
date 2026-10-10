import numpy as np
from numpy.polynomial.hermite_e import hermeval
from scipy.stats import norm
# alpha_k across final neurons ~ N(0, rho/(1-rho)) with rho=rho_15 ~ 0.916 (trajectory) ; sigma ~ 0.4
rho=0.916; v=rho/(1-rho); sig=0.40
a=np.linspace(-40,40,400001); wgt=norm.pdf(a,scale=np.sqrt(v)); da=a[1]-a[0]
def rms_He(j):
    c=np.zeros(j+1); c[j]=1
    f=hermeval(-a,c)*norm.pdf(a)
    return np.sqrt(np.sum(f**2*wgt)*da)
# standardised cumulant sizes measured at n=1024, L=16 (rms over neurons)
sk=0.116; k4c=0.066; k4i=0.026
terms={ 'kappa3 (m3)':(3,sk,1/6), 'kappa4 coherent (m4)':(4,k4c,1/24), 'kappa4 incoherent (m4)':(4,k4i,1/24),
        'm6=10 kappa3^2':(6,10*sk**2,1/720), 'm7=35 kappa3 kappa4 (total k4)':(7,35*sk*k4c,1/5040),
        'm7=35 kappa3 kappa4^inc':(7,35*sk*k4i,1/5040), 'm8=35 kappa4^2 (total)':(8,35*k4c**2,1/40320),
        'm9 ~ 280 kappa3^3':(9,280*sk**3,1/362880)}
print("rms He_{r-2}(-a)phi(a):", {r:round(rms_He(r-2),4) for r in range(3,10)})
for name,(r,m,c) in terms.items():
    print(f"{name:34s}: standardized size {m:.4f} -> rms mean shift sigma*m/r!*|He_(r-2) phi| = {sig*m*c*rms_He(r-2):.2e}")
print("target RMS 6e-5")
