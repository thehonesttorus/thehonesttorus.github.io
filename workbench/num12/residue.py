# Weights-only residue: the gain (scale-mixture) variance gamma, injected by each layer at tree level
# (coherent two-point fourth cumulant, averaged over pairs, O(n^2)) and transported critically (factor 1).
# Mean correction: m <- (1 - gamma/8) m_closure.
import numpy as np
from math import factorial
from closure import relu_coeffs, relu2_coeffs, relu_var
from edgeworth import relu_central
def closure_with_residue(Ws, K=14, terms=("k4", "F0", "k3mu", "Kpmu"), tau=1.0):
    n = Ws[0].shape[0]; out = []; gam = 0.0
    for l, W in enumerate(Ws):
        if l == 0:
            mu = np.zeros(n); S = W @ W.T
        else:
            mu = W @ m; S = W @ C @ W.T
            # fresh gain injected by this layer from the Gaussian source (prev)
            pmu, psig, pR, pm = prev
            q = mu**2 + np.diag(S)
            k3h, k4h = relu_central(pmu, psig)
            W2q = (W*W)/q[:, None]
            u = W2q.sum(0); v = (W2q**2).sum(0)
            t = W.T @ (mu/q)
            A = relu_coeffs(pmu, psig, K); B = relu2_coeffs(pmu, psig, K)
            c = B - 2*A[0]*A; c[0] = 0
            R0 = pR.copy(); np.fill_diagonal(R0, 0.0)
            Cov = np.zeros((n, n)); F = np.zeros((n, n)); Kp = np.zeros((n, n)); Rk = np.ones((n, n))
            for j in range(1, K+1):
                Rk = Rk*R0; f = Rk/factorial(j)
                Cov += np.outer(A[j], A[j])*f; F += np.outer(c[j], c[j])*f; Kp += np.outer(c[j], A[j])*f
            F0 = F - 2*Cov*Cov
            tot = 0.0
            if "k4" in terms: tot += k4h @ (u*u - v)
            if "F0" in terms: tot += u @ F0 @ u
            if "k3mu" in terms: tot += 4*np.sum(k3h*u*t)
            if "Kpmu" in terms: tot += 4*(u @ Kp @ t)
            gam = tau*gam + tot/(n*(n-1))
        sig = np.sqrt(np.diag(S)); R = S/np.outer(sig, sig)
        A = relu_coeffs(mu, sig, K)
        m = A[0]; C = np.zeros((n, n)); Rk = np.ones((n, n))
        for k in range(1, K+1):
            Rk = Rk*R; C += np.outer(A[k], A[k])*Rk/factorial(k)
        np.fill_diagonal(C, relu_var(mu, sig))
        prev = (mu, sig, R, m)
        out.append(dict(mu=mu, sig=sig, m=m*(1-gam/8), m_closure=m.copy(), gamma=gam))
    return out
