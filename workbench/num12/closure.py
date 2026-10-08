# Gaussian closure with the actual weights: propagate mean vector and full covariance.
import numpy as np
from scipy.special import ndtr
from math import factorial
SQ2PI = np.sqrt(2*np.pi)
def phi(a): return np.exp(-0.5*a*a)/SQ2PI
def herm(a, K):
    # probabilists' Hermite He_0..He_K at array a
    H = [np.ones_like(a), a.copy()]
    for k in range(1, K):
        H.append(a*H[k] - k*H[k-1])
    return H
def relu_coeffs(mu, sig, K):
    """a_j = E[g(z) He_j(u)], g = ReLU, z = mu + sig u; j = 0..K."""
    a = mu/sig; P = ndtr(a); f = phi(a); H = herm(a, K)
    A = np.zeros((K+1,) + a.shape)
    A[0] = sig*(a*P + f); A[1] = sig*P
    for j in range(2, K+1):
        A[j] = sig*((-1)**j)*H[j-2]*f
    return A
def relu2_coeffs(mu, sig, K):
    """b_j = E[g(z)^2 He_j(u)]."""
    a = mu/sig; P = ndtr(a); f = phi(a); H = herm(a, K)
    B = np.zeros((K+1,) + a.shape)
    B[0] = sig**2*((1+a*a)*P + a*f); B[1] = sig**2*2*(a*P + f); B[2] = sig**2*2*P
    for j in range(3, K+1):
        B[j] = sig**2*2*((-1)**(j-3))*H[j-3]*f
    return B
def relu_var(mu, sig):
    a = mu/sig; P = ndtr(a); f = phi(a)
    M = sig*(a*P + f)
    return sig**2*((1+a*a)*P + a*f) - M*M
def closure(Ws, K=14, keep=False):
    n = Ws[0].shape[0]
    mu = np.zeros(n); S = Ws[0] @ Ws[0].T
    out = []
    for l, W in enumerate(Ws):
        if l > 0:
            mu = W @ m; S = W @ C @ W.T
        sig = np.sqrt(np.diag(S)); R = S/np.outer(sig, sig)
        A = relu_coeffs(mu, sig, K)
        m = A[0]
        C = np.zeros((n, n)); Rk = np.ones((n, n))
        for k in range(1, K+1):
            Rk = Rk*R
            C += np.outer(A[k], A[k])*Rk/factorial(k)
        np.fill_diagonal(C, relu_var(mu, sig))
        out.append(dict(mu=mu, sig=sig, R=R, m=m, C=C) if keep else dict(mu=mu, sig=sig, m=m))
    return out
