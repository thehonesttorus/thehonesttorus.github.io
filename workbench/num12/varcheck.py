import numpy as np, sys
from closure import closure, relu_var
from corrected import corrected_closure, shift_coeffs
from edgeworth import cumulants_next
n, L, s = 256, 4, 0
Ws = list(np.load(f"W_n{n}_L{L}_s{s}.npy")); tr = np.load(f"truth_n{n}_L{L}_s{s}.npz")
m, v = tr["m"], tr["v"]
out = closure(Ws, keep=True)
o1, o2 = out[0], out[1]
k3, k4, _ = cumulants_next(Ws[1], o1["mu"], o1["sig"], o1["R"])
mu, sig = o2["mu"], o2["sig"]
a3,a4,a6,b3,b4,b6 = shift_coeffs(mu, sig)
s3, s4 = k3/(6*sig**3), k4/(24*sig**4)
m0 = o2["m"]; v0 = relu_var(mu, sig)
dm = s3*a3+s4*a4+s3*s3*0.5*a6; d2 = s3*b3+s4*b4+s3*s3*0.5*b6
v1 = v0+m0*m0+d2-(m0+dm)**2
print("var err closure %.3e corrected %.3e  (rel var noise ~ %.1e)" % (np.sqrt(np.mean((v0-v[1])**2)), np.sqrt(np.mean((v1-v[1])**2)), np.sqrt(2/tr['T'])*v[1].mean()))
print("mean var err closure %+.3e corrected %+.3e" % ((v0-v[1]).mean(), (v1-v[1]).mean()))
# check E g^2 Hermite coefficient numerically
from scipy.integrate import quad
from scipy.special import eval_hermitenorm as He
from closure import phi
mu_,sg_=0.3,1.2; a=mu_/sg_
for j,b in [(3,b3),(4,b4),(6,b6)]:
    num = quad(lambda u:(mu_+sg_*u)**2*He(j,u)*phi(u),-a,40)[0]
    _,_,_,B3,B4,B6 = shift_coeffs(np.array([mu_]),np.array([sg_]))
    print(j, num, {3:B3,4:B4,6:B6}[j])
