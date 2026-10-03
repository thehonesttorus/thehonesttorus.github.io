# Budget version: tree calculus restricted to single-neuron (diagonal) sources, which cost one n^3 product
# per two-point term. Sources: kappa3(h_a), kappa4(h_a). Optional pair terms PP4 (annealed) and response r<=rmax.
import numpy as np
from math import factorial
from scipy.special import ndtr
from closure import relu_coeffs, relu2_coeffs, relu_var
from edgeworth import relu_central
from corrected import shift_coeffs
from corrected3 import mehler_shift
def diag_tlp(Ws, K=12, rmax=1, k21=True, k22=True, counter=None):
    n = Ws[0].shape[0]; out = []; hist = []
    mm = lambda X, Y: (counter.__setitem__(0, counter[0]+1) if counter is not None else None) or X @ Y
    for l, W in enumerate(Ws):
        if l == 0:
            mu = np.zeros(n); S = mm(W, W.T)
        else:
            mu = W @ m; S = mm(mm(W, C), W.T)
        sig = np.sqrt(np.diag(S)); R = S/np.outer(sig, sig)
        A = relu_coeffs(mu, sig, K+6)
        m = A[0].copy(); C = np.zeros((n, n)); Rk = np.ones((n, n))
        for k in range(1, K+1):
            Rk = Rk*R; C += np.outer(A[k], A[k])*Rk/factorial(k)
        v = relu_var(mu, sig)
        if l > 0:
            k3 = np.zeros(n); k4 = np.zeros(n); K21 = np.zeros((n, n)); K22 = np.zeros((n, n)); B = W
            for r in range(0, min(l, rmax+1)):
                k3h, k4h, P = hist[l-1-r]
                B2 = B*B
                k3 += (B2*B) @ k3h; k4 += (B2*B2) @ k4h
                if k21: K21 += mm(B2*k3h, B.T)
                if k22: K22 += mm(B2*k4h, B2.T)
                if l-2-r >= 0 and r < rmax: B = mm(B*P, Ws[l-1-r])
            a3, a4, a6, b3, b4, b6 = shift_coeffs(mu, sig)
            s3, s4 = k3/(6*sig**3), k4/(24*sig**4)
            dm = s3*a3 + s4*a4 + 0.5*s3*s3*a6; d2 = s3*b3 + s4*b4 + 0.5*s3*s3*b6
            E2 = v + m*m; m = m + dm; v = E2 + d2 - m*m
            dC = 0
            if k21: dC = dC + (K21/2)*mehler_shift(A, sig, 2, 1, R) + (K21.T/2)*mehler_shift(A, sig, 1, 2, R)
            if k22: dC = dC + (K22/4)*mehler_shift(A, sig, 2, 2, R)
            dC = dC + (k3[:, None]/6)*mehler_shift(A, sig, 3, 0, R, pmin=1) + (k3[None, :]/6)*mehler_shift(A, sig, 0, 3, R, pmin=1)
            dC = dC + (k4[:, None]/24)*mehler_shift(A, sig, 4, 0, R, pmin=1) + (k4[None, :]/24)*mehler_shift(A, sig, 0, 4, R, pmin=1)
            C = C + dC
        np.fill_diagonal(C, v)
        k3h, k4h = relu_central(mu, sig)
        hist.append((k3h, k4h, ndtr(mu/sig)))
        out.append(dict(mu=mu, sig=sig, m=m))
    return out
