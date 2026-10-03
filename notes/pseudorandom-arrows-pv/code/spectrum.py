# Taylor coefficients b_r^{(D)} of the D-fold iterate of the ReLU correlation map
# rho(c) = c/2 + (sqrt(1-c^2) + c*arcsin c)/pi   (He init, bias-free ReLU)
import numpy as np
from math import comb, pi
R = 600  # truncation degree
# series of sqrt(1-c^2) + c*arcsin(c) = 1 + sum_{k>=1} c^{2k} * [(2k-3)!!]^2 / (2k)!  ... compute directly
a = np.zeros(R+1)
# sqrt(1-x) = sum_k binom(1/2,k)(-x)^k ; arcsin c = sum_k (2k)!/(4^k (k!)^2 (2k+1)) c^{2k+1}
from fractions import Fraction
def binom_half(k):
    r = Fraction(1)
    for i in range(k):
        r *= Fraction(1,2) - i
        r /= i+1
    return r
for k in range(0, R//2+1):
    s = float(binom_half(k)) * (-1)**k           # sqrt(1-c^2) coefficient of c^{2k}
    if 2*k <= R: a[2*k] += s
    if k >= 1:                                    # c*arcsin c coefficient of c^{2k}
        j = k-1
        coef = comb(2*j, j) / (4**j * (2*j+1))
        a[2*k] += coef
a /= pi
a[1] += 0.5
print("rho coefficients r=0..6:", np.round(a[:7], 6), " sum(trunc) =", a.sum())
def compose(p, q):
    # p(q(c)) truncated, Horner
    out = np.zeros(R+1); out[0] = p[-1]
    for coef in p[::-1][1:]:
        out = np.convolve(out, q)[:R+1]
        out[0] += coef
    return out
b = a.copy()
spec = {1: b.copy()}
for D in range(2, 17):
    b = compose(a, b)   # rho(rho^{D-1})
    spec[D] = b.copy()
print("D  b0      b1      b2      b3      b4      sum_{r>=5}  (trunc mass)")
for D in [1,2,3,4,6,8,12,15,16]:
    s = spec[D]
    print(f"{D:2d} {s[0]:.4f}  {s[1]:.4f}  {s[2]:.4f}  {s[3]:.4f}  {s[4]:.4f}  {s[5:].sum():.4f}     {s.sum():.4f}")
np.save("spec.npy", np.array([spec[D] for D in range(1,17)]))
