# Closure + tree-level cumulants + linear-response transport of earlier non-Gaussianity.
import numpy as np
from math import factorial
from scipy.special import ndtr
from closure import relu_coeffs, relu_var, phi
from edgeworth import cumulants_next
from corrected import shift_coeffs
def corrected_closure2(Ws, K=14, rmax=99, var_corr=True, J=4):
    n = Ws[0].shape[0]; out = []; hist = []   # hist[j] = (mu, sig, R, P) of pre-activation at layer j
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
            k3 = np.zeros(n); k4 = np.zeros(n); B = W
            for r in range(0, min(l, rmax+1)):
                pm, ps, pR, pP = hist[l-1-r]
                t3, t4, _ = cumulants_next(B, pm, ps, pR, J=J)
                k3 += t3; k4 += t4
                if l-2-r >= 0:
                    B = (B*pP) @ Ws[l-1-r]
            a3, a4, a6, b3, b4, b6 = shift_coeffs(mu, sig)
            s3, s4 = k3/(6*sig**3), k4/(24*sig**4)
            dm = s3*a3 + s4*a4 + s3*s3*0.5*a6
            d2 = s3*b3 + s4*b4 + s3*s3*0.5*b6
            E2 = v + m*m
            m = m + dm
            if var_corr: v = E2 + d2 - m*m
        np.fill_diagonal(C, v)
        hist.append((mu, sig, R, ndtr(mu/sig)))
        out.append(dict(mu=mu, sig=sig, m=m))
    return out
