# Two-point tree-level cumulants K21_kl = k(y_k,y_k,y_l), K22_kl = k(y_k,y_k,y_l,y_l) for y = B g(x),
# x Gaussian (closure data mu, sig, R at the source layer). Off-diagonal entries only are meaningful.
import numpy as np
from math import factorial
from closure import relu_coeffs, relu2_coeffs
from edgeworth import relu_central
def twopoint(B, mu, sig, R, J=3):
    n = len(mu); K = 2*J+2
    A = relu_coeffs(mu, sig, K); Bc = relu2_coeffs(mu, sig, K)
    c = Bc - 2*A[0]*A; c[0] = 0
    R0 = R.copy(); np.fill_diagonal(R0, 0.0)
    Rp = [np.ones_like(R0)]
    for p in range(1, K+1): Rp.append(Rp[-1]*R0)
    B2 = B*B
    k3h, k4h = relu_central(mu, sig)
    # kappa(h_a,h_a,h_c) and centred-square Mehler matrices
    Kp = np.zeros_like(R0); Cov = np.zeros_like(R0); F = np.zeros_like(R0)
    for j in range(1, K+1):
        fj = Rp[j]/factorial(j)
        Kp += np.outer(c[j], A[j])*fj; Cov += np.outer(A[j], A[j])*fj; F += np.outer(c[j], c[j])*fj
    F0 = F - 2*Cov*Cov; np.fill_diagonal(F0, 0.0)
    K21 = (B2*k3h) @ B.T + B2 @ Kp @ B.T + 2*((B*(B @ Kp.T)) @ B.T)
    U = {p: (B*A[p]) @ (Rp[p]/factorial(p)) for p in range(1, J+1)}
    for p in range(1, J+1):
        for q in range(1, J+1):
            fpq = factorial(p)*factorial(q)
            # centre at a (and b): 2 * sum_a B_ka a_{a,p+q} U^p[k,a] U^q[l,a]  minus b=c coincidence
            cA = (B*A[p+q]*U[p]) @ U[q].T
            Z = B*(((B*A[p+q]) @ Rp[p+q])*(A[p]*A[q])/fpq)
            cA -= Z @ B.T
            # centre at c: sum_c B_lc a_{c,p+q} U^p[k,c] U^q[k,c]  minus a=b coincidence
            V = (B2*(A[p]*A[q])) @ (Rp[p+q]/fpq)
            cC = (U[p]*U[q] - V) @ (B*A[p+q]).T
            K21 += 2*cA + cC
    e = c[2] - 2*A[1]**2
    V11 = (B2*A[1]**2) @ Rp[2]
    G = U[1]**2 - V11
    K22 = (B2*k4h) @ B2.T + B2 @ F0 @ B2.T + (B2*e) @ G.T + G @ (B2*e).T
    return K21, K22
