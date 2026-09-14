# The Mac report's minimal example, checked analytically end to end.
#   X ~ N(0,1),  rho(rho(X)) = rho(X),  so the exact final mean response to ANY perturbation is the
#   response of a single ReLU.  For the covariance-matched cubic direction R(v,Sigma) with v=1,
#   Sigma=1 the tensor is R(1,1) = 3, i.e. a kappa_3 perturbation of size 3.
# First-order cumulant response:  d E[g(X)] = (kappa_3/6) E[g'''(X)].
import numpy as np
from scipy.special import ndtr
phi=lambda t: np.exp(-0.5*t*t)/np.sqrt(2*np.pi)
k3=3.0                                   # R(1,1)
# (rho)'''  = delta'   ->  E[delta'(X)] = -p'(0) = 0 at a CENTRED Gaussian
d_mean_exact = (k3/6.0)*0.0
# (rho^2)''' = 2 delta  ->  E[2 delta(X)] = 2 phi(0)
d_var_after1 = (k3/6.0)*2*phi(0.0)       # d E[rho^2]; d(mean)=0 so dVar = dE[rho^2]
print("exact final mean response                      : %.8f   (zero: rho o rho = rho)" % d_mean_exact)
print("variance response after one ReLU               : %.8f" % d_var_after1)
print("  report says                                  : 0.39894228")
# Now the closure: replace the law of rho(X) by a Gaussian with its mean and PERTURBED variance,
# then apply the Gaussian-ReLU mean formula Psi(m,s) = m Phi(m/s) + s phi(m/s).  dPsi/ds = phi(m/s).
m1=phi(0.0); v1=0.5-1/(2*np.pi); s1=np.sqrt(v1)
t1=m1/s1
d_s = d_var_after1/(2*s1)                # ds = dvar/(2s)
d_mean_closure = phi(t1)*d_s             # dm = 0 at the centred reference
print()
print("closure intermediate: m=%.8f  var=%.8f  s=%.8f  t=%.8f" % (m1,v1,s1,t1))
print("closure-induced final mean response            : %.8f" % d_mean_closure)
print("  report says                                  : 0.10792360")
print()
print("ratio spurious/true variance-response gain     : %.6f" % (phi(t1)/(2*s1)))
# finite-difference confirmation of the closure step, no derivatives
eps=1e-5
def Psi(m,s): 
    t=m/s; return m*ndtr(t)+s*phi(t)
fd=(Psi(m1,np.sqrt(v1+eps*d_var_after1))-Psi(m1,np.sqrt(v1-eps*d_var_after1)))/(2*eps)
print("finite-difference check of the same quantity   : %.8f" % fd)
