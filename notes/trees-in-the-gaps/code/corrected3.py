# Weights-only tree-level cumulant propagation: closure + one-point (mean/var) and two-point (covariance)
# Edgeworth corrections, with sources transported linearly from up to rmax earlier layers.
import numpy as np
from math import factorial
from scipy.special import ndtr
from closure import relu_coeffs, relu_var
from edgeworth import cumulants_next
from corrected import shift_coeffs
from twopoint import twopoint
from blob import g1_kappa3
from closure import phi
def mehler_shift(A, sig, i, j, R, pmin=0, P=8):
    """M_kl = sum_{p>=pmin} sig_k^-i sig_l^-j a_{k,i+p} a_{l,j+p} R^p/p!"""
    M = np.zeros_like(R); Rp = np.ones_like(R)
    for p in range(0, P+1):
        if p > 0: Rp = Rp*R
        if p >= pmin:
            M += np.outer(A[i+p]/sig**i, A[j+p]/sig**j)*Rp/factorial(p)
    return M
def tree_closure(Ws, K=14, rmax=2, J=3, two_point=True, one_point=True, var_corr=True, gap_blob=False):
    n = Ws[0].shape[0]; out = []; hist = []
    for l, W in enumerate(Ws):
        if l == 0:
            mu = np.zeros(n); S = W @ W.T
        else:
            mu = W @ m; S = W @ C @ W.T
        sig = np.sqrt(np.diag(S)); R = S/np.outer(sig, sig)
        A = relu_coeffs(mu, sig, K+8)
        m = A[0].copy()
        C = np.zeros((n, n)); Rk = np.ones((n, n))
        for k in range(1, K+1):
            Rk = Rk*R; C += np.outer(A[k], A[k])*Rk/factorial(k)
        v = relu_var(mu, sig)
        if l > 0:
            k3 = np.zeros(n); k4 = np.zeros(n); K21 = np.zeros((n, n)); K22 = np.zeros((n, n)); B = W
            for r in range(0, min(l, rmax+1)):
                pm, ps, pR, pP = hist[l-1-r]
                t3, t4, _ = cumulants_next(B, pm, ps, pR, J=J); k3 += t3; k4 += t4
                if two_point:
                    s21, s22 = twopoint(B, pm, ps, pR, J=J); K21 += s21; K22 += s22
                if l-2-r >= 0: B = (B*pP) @ Ws[l-1-r]
            if gap_blob and l >= 2:
                pm1, ps1, pR1, pP1 = hist[l-1]; pm2, ps2, pR2, _ = hist[l-2]
                k3 = k3 + g1_kappa3(W, Ws[l-1], pP1, phi(pm1/ps1)/ps1, pm2, ps2, pR2)
            a3, a4, a6, b3, b4, b6 = shift_coeffs(mu, sig)
            s3, s4 = k3/(6*sig**3), k4/(24*sig**4)
            dm = s3*a3 + s4*a4 + s3*s3*0.5*a6
            d2 = s3*b3 + s4*b4 + s3*s3*0.5*b6
            E2 = v + m*m
            if one_point:
                m = m + dm
                if var_corr: v = E2 + d2 - m*m
            if two_point:
                dC = (K21/2)*mehler_shift(A, sig, 2, 1, R) + (K21.T/2)*mehler_shift(A, sig, 1, 2, R) \
                     + (K22/4)*mehler_shift(A, sig, 2, 2, R)
                dC += (k3[:, None]/6)*mehler_shift(A, sig, 3, 0, R, pmin=1) + (k3[None, :]/6)*mehler_shift(A, sig, 0, 3, R, pmin=1)
                dC += (k4[:, None]/24)*mehler_shift(A, sig, 4, 0, R, pmin=1) + (k4[None, :]/24)*mehler_shift(A, sig, 0, 4, R, pmin=1)
                C += dC
        np.fill_diagonal(C, v)
        hist.append((mu, sig, R, ndtr(mu/sig)))
        out.append(dict(mu=mu, sig=sig, m=m))
    return out
