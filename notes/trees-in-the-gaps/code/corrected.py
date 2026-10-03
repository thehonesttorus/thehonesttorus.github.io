# Closure + tree-level Edgeworth corrections to means and variances, layer by layer.
import numpy as np
from math import factorial
from closure import relu_coeffs, relu2_coeffs, relu_var, phi
from edgeworth import cumulants_next
def shift_coeffs(mu, sig):
    a = mu/sig; f = phi(a)
    a3 = -sig*a*f; a4 = sig*(a*a-1)*f; a6 = sig*(a**4-6*a*a+3)*f
    # b_j = E[g^2 He_j] = 2 sig^2 (-1)^(j-3) He_{j-3}(a) phi(a) for j>=3
    b3 = 2*sig**2*f; b4 = -2*sig**2*a*f; b6 = -2*sig**2*(a**3-3*a)*f
    return a3, a4, a6, b3, b4, b6
def corrected_closure(Ws, K=14, var_corr=True, mean_corr=True):
    n = Ws[0].shape[0]; out = []
    prev = None
    for l, W in enumerate(Ws):
        if l == 0:
            mu = np.zeros(n); S = W @ W.T
        else:
            mu = W @ m; S = W @ C @ W.T
        sig = np.sqrt(np.diag(S)); R = S/np.outer(sig, sig)
        A = relu_coeffs(mu, sig, K)
        m = A[0].copy()
        C = np.zeros((n, n)); Rk = np.ones((n, n))
        for k in range(1, K+1):
            Rk = Rk*R; C += np.outer(A[k], A[k])*Rk/factorial(k)
        v = relu_var(mu, sig)
        if l > 0:
            k3, k4, _ = cumulants_next(W, *prev)
            a3, a4, a6, b3, b4, b6 = shift_coeffs(mu, sig)
            s3, s4 = k3/(6*sig**3), k4/(24*sig**4)
            dm = s3*a3 + s4*a4 + s3*s3*0.5*a6   # k3^2/72 = (k3/6)^2/2
            d2 = s3*b3 + s4*b4 + s3*s3*0.5*b6
            E2 = v + m*m
            if mean_corr: m = m + dm
            if var_corr: v = E2 + d2 - m*m
        np.fill_diagonal(C, v)
        prev = (mu, sig, R)
        out.append(dict(mu=mu, sig=sig, m=m))
    return out
