# TLP-lite: closure + the O(n^2)-per-layer part of the tree calculus (quenched diagonal D3, coherent
# fourth-cumulant terms annealed over the source weights, coherent two-point K22, annealed response).
import numpy as np
from math import factorial
from scipy.special import ndtr
from closure import relu_coeffs, relu2_coeffs, relu_var, phi
from edgeworth import relu_central
from corrected import shift_coeffs
def lite(Ws, K=12, rmax=4, use=("D3", "D4", "PP4", "K22", "resp")):
    n = Ws[0].shape[0]; out = []; src = []   # src[j] = (sum kappa4, Fbar, chibar) of activation layer j
    m = C = None
    for l, W in enumerate(Ws):
        if l == 0:
            mu = np.zeros(n); S = W @ W.T
        else:
            mu = W @ m; S = W @ C @ W.T
        sig = np.sqrt(np.diag(S)); R = S/np.outer(sig, sig)
        A = relu_coeffs(mu, sig, K+6); Bc = relu2_coeffs(mu, sig, K)
        R0 = R.copy(); np.fill_diagonal(R0, 0.0)
        Cm = np.zeros((n, n)); F = np.zeros((n, n)); Rk = np.ones((n, n))
        c = Bc - 2*A[0]*A[:K+1]; c[0] = 0
        for k in range(1, K+1):
            Rk = Rk*R0; f = Rk/factorial(k)
            Cm += np.outer(A[k], A[k])*f; F += np.outer(c[k], c[k])*f
        m = A[0].copy(); v = relu_var(mu, sig)
        if l > 0:
            k3h, k4h, s4, Fbar, chib = prev
            rn = (W*W).sum(1)                          # row norms ||W_i||^2
            k3 = (W**3) @ k3h if "D3" in use else 0*m
            k4 = (W**4) @ k4h if "D4" in use else 0*m
            if "PP4" in use: k4 = k4 + 3*(rn/n)**2*Fbar
            Kbar = (s4 + Fbar)/n**2                    # annealed K22 per unit row norm^2
            if "resp" in use:
                b2 = rn.copy()
                for r in range(1, min(l, rmax+1)):
                    b2 = b2*src[l-r][2]                # ||B^{(r)}_i||^2 ~ chibar * ||B^{(r-1)}_i||^2
                    s4r, Fr = src[l-1-r][0], src[l-1-r][1]
                    k4 = k4 + 3*(b2/n)**2*(s4r + Fr)
                    Kbar = Kbar + (b2.mean()/n)**2*(s4r + Fr)/rn.mean()**2*1.0
            a3, a4, a6, b3, b4, b6 = shift_coeffs(mu, sig)
            s3_, s4_ = k3/(6*sig**3), k4/(24*sig**4)
            dm = s3_*a3 + s4_*a4 + 0.5*s3_**2*a6
            d2 = s3_*b3 + s4_*b4 + 0.5*s3_**2*b6
            E2 = v + m*m; m = m + dm; v = E2 + d2 - m*m
            if "K22" in use:
                kk = Kbar*np.outer(rn, rn)             # K22_kl ~ ||W_k||^2 ||W_l||^2 (sum k4 + Fbar)/n^2
                M = np.zeros((n, n)); Rp = np.ones((n, n))
                for p in range(0, 6):
                    if p > 0: Rp = Rp*R
                    M += np.outer(A[2+p]/sig**2, A[2+p]/sig**2)*Rp/factorial(p)
                Cm += kk/4*M
        C = Cm; np.fill_diagonal(C, v)
        k3h, k4h = relu_central(mu, sig)
        F0 = F - 2*Cm*Cm; np.fill_diagonal(F0, 0.0)
        P = ndtr(mu/sig)
        prev = (k3h, k4h, k4h.sum(), F0.sum(), 2*np.mean(P**2))
        src.append(prev[2:4] + (prev[4],))
        out.append(dict(mu=mu, sig=sig, m=m))
    return out
