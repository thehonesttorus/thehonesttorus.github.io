# One-blob diagrams: an earlier layer's 4-point blob attached through a gap vertex.
import numpy as np
from math import factorial
from closure import relu_coeffs, relu2_coeffs
from edgeworth import relu_central
def k22_cross(U, V, mu, sig, R, K=8):
    """M[a,i] = kappa(u_a.h, u_a.h, v_i.h, v_i.h), h = relu(x), x ~ N(mu, Sigma) (tree level, coherent terms)."""
    A = relu_coeffs(mu, sig, K); Bc = relu2_coeffs(mu, sig, K)
    c = Bc - 2*A[0]*A; c[0] = 0
    R0 = R.copy(); np.fill_diagonal(R0, 0.0)
    Rp = [np.ones_like(R0)]
    for p in range(1, K+1): Rp.append(Rp[-1]*R0)
    k3h, k4h = relu_central(mu, sig)
    Cov = np.zeros_like(R0); F = np.zeros_like(R0)
    for j in range(1, K+1):
        fj = Rp[j]/factorial(j); Cov += np.outer(A[j], A[j])*fj; F += np.outer(c[j], c[j])*fj
    F0 = F - 2*Cov*Cov; np.fill_diagonal(F0, 0.0)
    U2, V2 = U*U, V*V
    e = c[2] - 2*A[1]**2
    def G(X):
        X2 = X*X; Ux = (X*A[1]) @ R0
        return Ux**2 - (X2*A[1]**2) @ Rp[2]
    M = (U2*k4h) @ V2.T + U2 @ F0 @ V2.T + (U2*e) @ G(V).T + G(U) @ (V2*e).T
    # cross pairing {a: k,l}{i: k,l}: 2 sum_{k!=l} u_k u_l v_k v_l F0_kl  (needs per-pair; use (U*V) trick)
    # sum_{k,l} (u_k v_k)(u_l v_l) F0_kl = rowpair bilinear -> compute via einsum on small n only
    return M
def g1_kappa3(Wnext, Wcur, Pcur, pcur, mu_prev, sig_prev, R_prev):
    """(3/2) sum_a Wnext_ia p_a kappa(y_a,y_a, beta_i.h, beta_i.h), y = Wcur h, beta_i = Wcur^T (Wnext_i * Pcur)."""
    Beta = (Wnext*Pcur) @ Wcur
    M = k22_cross(Wcur, Beta, mu_prev, sig_prev, R_prev)
    return 1.5*np.sum(Wnext*pcur*M.T, 1)
