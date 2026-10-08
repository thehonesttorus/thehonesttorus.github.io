# Tree-level cumulants of z = W g(y), y Gaussian (closure data), and Edgeworth mean corrections.
import numpy as np
from math import factorial
from scipy.special import ndtr
from closure import relu_coeffs, relu2_coeffs, relu_var, phi
def relu_central(mu, sig):
    a = mu/sig; P = ndtr(a); f = phi(a)
    E1 = sig*(a*P+f); E2 = sig**2*((1+a*a)*P+a*f)
    E3 = sig**3*((a**3+3*a)*P+(a*a+2)*f); E4 = sig**4*((a**4+6*a*a+3)*P+(a**3+5*a)*f)
    m = E1; v = E2-m*m
    k3 = E3-3*m*E2+2*m**3
    k4 = E4-4*m*E3+6*m*m*E2-3*m**4-3*v*v
    return k3, k4
def cumulants_next(W, mu, sig, R, J=4, terms=("D3","P3","T3","D4","PP4","PD4")):
    """kappa_3, kappa_4 of z_i = sum_k W_ik g(y_k), y ~ N(mu, Sigma) with correlations R."""
    n = len(mu); K = 2*J+2
    A = relu_coeffs(mu, sig, K); B = relu2_coeffs(mu, sig, K)
    c = B - 2*A[0]*A; c[0] = 0
    R0 = R.copy(); np.fill_diagonal(R0, 0.0)
    Rp = [np.ones_like(R0)]
    for p in range(1, K+1): Rp.append(Rp[-1]*R0)
    W2 = W*W
    k3h, k4h = relu_central(mu, sig)
    out = {}
    out["D3"] = (W2*W) @ k3h
    P3 = np.zeros(n)
    for j in range(1, K+1):
        X = (W2*c[j]) @ (Rp[j]/factorial(j)); P3 += 3*np.sum(X*(W*A[j]), 1)
    out["P3"] = P3
    U = {p: (W*A[p]) @ (Rp[p]/factorial(p)) for p in range(1, J+1)}
    T3 = np.zeros(n)
    for p in range(1, J+1):
        for q in range(1, J+1):
            V = (W2*(A[p]*A[q])) @ (Rp[p+q]/(factorial(p)*factorial(q)))
            T3 += 3*np.sum(W*A[p+q]*(U[p]*U[q]-V), 1)
    out["T3"] = T3
    out["D4"] = (W2*W2) @ k4h
    Cov = np.zeros_like(R0); F = np.zeros_like(R0)
    for j in range(1, K+1):
        Cov += np.outer(A[j], A[j])*Rp[j]/factorial(j); F += np.outer(c[j], c[j])*Rp[j]/factorial(j)
    F0 = F-2*Cov*Cov; np.fill_diagonal(F0, 0.0)
    out["PP4"] = 3*np.sum((W2 @ F0)*W2, 1)
    e = c[2]-2*A[1]**2
    V11 = (W2*A[1]**2) @ Rp[2]
    Q = (W2*c[1]) @ Rp[1]
    out["PD4"] = 6*np.sum(W2*e*(U[1]**2-V11), 1) + 12*np.sum(W*A[2]*U[1]*Q, 1)
    k3 = sum(out[t] for t in terms if t.endswith("3")); k4 = sum(out[t] for t in terms if t.endswith("4"))
    return k3, k4, out
def edgeworth_shift(mu, sig, k3, k4):
    a = mu/sig; f = phi(a)
    a3 = -sig*a*f; a4 = sig*(a*a-1)*f; a6 = sig*(a**4-6*a*a+3)*f
    return k3/(6*sig**3)*a3 + k4/(24*sig**4)*a4 + k3**2/(72*sig**6)*a6
