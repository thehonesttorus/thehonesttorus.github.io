# Gain-aware closure + fresh kappa_3 trees, flopscope, official orientation (z = V^T h).  Fewer ops than gac_fs.py:
# the relu^2 Hermite coefficients are B_j = 2 sigma A_{j-1} (j >= 1), so the previous layer's A is reused everywhere.
from math import lgamma, log, exp, factorial
import flopscope as flops
import flopscope.numpy as fnp
F32 = fnp.float32
def EG(g):
    if g <= 0: return 1.0
    k = 1.0/g; return exp(lgamma(k + 0.5) - lgamma(k) - 0.5*log(k))
def predict(weights, trees=True, K=14):
    n = weights[0].shape[0]; rows = []; gam = 0.0; g1 = 1.0; zero = fnp.zeros((), dtype=F32)
    rf = [1.0/factorial(j) for j in range(K+1)]; srf = [x**0.5 for x in rf]
    for l, V in enumerate(weights):
        if l == 0:
            S = fnp.matmul(V.T, V); sv = fnp.diag(S); mu = fnp.zeros(n, dtype=F32)
        else:
            w = V.T @ M; Qc = fnp.einsum("ij,ia,jb->ab", Csym, V, V); q = fnp.diag(Qc) + w*w
            iq = 1.0/q; V2 = V*V
            u = V2 @ iq; v = (V2*V2) @ (iq*iq); t = V @ (w*iq)
            # previous conditional state: A (list), sig, a, P, f, R0 (zero diagonal), C (closure covariance)
            c = [None] + [2*(psig*pA[j-1] - pA[0]*pA[j]) for j in range(1, K+1)]
            Cv = Cprev.copy(); fnp.fill_diagonal(Cv, zero)
            tot = fnp.sum(k4h*(u*u - v)) - 2*fnp.dot(u, (Cv*Cv) @ u) + 4*fnp.sum(k3h*u*t)
            Rj = R0; X = fnp.stack([u*c[1], t*pA[1]], axis=1)
            for j in range(1, K+1):
                if j > 1:
                    Rj = Rj*R0; X = fnp.stack([u*c[j], t*pA[j]], axis=1)
                Y = Rj @ X
                tot = tot + rf[j]*fnp.dot(X[:, 0], Y[:, 0] + 4*Y[:, 1])
                if j == 2: R02 = Rj
            g1p = g1; gam = gam + float(tot)/(n*(n - 1)); g1 = EG(gam); r = g1p/g1
            mu = r*w; S = Qc + (1.0 - r*r)*fnp.outer(w, w); sv = fnp.diag(S)
        sig = fnp.sqrt(fnp.maximum(sv, 1e-12)); R = S/fnp.outer(sig, sig)
        a = mu/sig; P = flops.stats.norm.cdf(a).astype(F32); f = flops.stats.norm.pdf(a).astype(F32)
        sf = sig*f; A = [sig*(a*P + f), sig*P, sf]; Hm, H = fnp.ones_like(a), a     # He_{j-2}: start He_0, He_1
        for j in range(3, K+1):
            A.append(((-1)**j)*sf*H); Hm, H = H, a*H - (j-2)*Hm
        C = fnp.outer(A[1], A[1])*R; Rk = R
        for k in range(2, K+1):
            Rk = Rk*R; s = A[k]*srf[k]; C = C + fnp.outer(s, s)*Rk
        M = A[0]; E2 = sig*sig*((1 + a*a)*P + a*f); fnp.fill_diagonal(C, E2 - M*M)
        if trees and l > 0:
            k3 = k3h @ (V2*V)
            Y1 = R0 @ (V2*c[1][:, None]); Y2 = R02 @ (V2*c[2][:, None])
            Ut = R0 @ (V*pA[1][:, None]); Vt = R02 @ (V2*(pA[1]*pA[1])[:, None])
            k3 = k3 + 3*fnp.sum(Y1*(V*pA[1][:, None]) + 0.5*Y2*(V*pA[2][:, None]) + (V*pA[2][:, None])*(Ut*Ut - Vt), axis=0)
            s3 = k3/(6*sig**3); a2 = a*a
            M = M + s3*(-sig*a*f) + 0.5*s3*s3*sig*(a2*a2 - 6*a2 + 3)*f
        rows.append(g1*M)
        # this layer becomes the source of the next: central moments of relu(y), y ~ N(mu, sig^2)
        E1 = A[0]; s2 = sig*sig; a2 = a*a
        E3 = s2*sig*((a2*a + 3*a)*P + (a2 + 2)*f); E4 = s2*s2*((a2*a2 + 6*a2 + 3)*P + (a2*a + 5*a)*f); vv = E2 - E1*E1
        k3h = E3 - 3*E1*E2 + 2*E1**3; k4h = E4 - 4*E1*E3 + 6*E1*E1*E2 - 3*E1**4 - 3*vv*vv
        pA, psig, Cprev = A, sig, C
        R0 = R.copy(); fnp.fill_diagonal(R0, zero)
        Csym = flops.as_symmetric(C, symmetry=(0, 1))
    return fnp.stack(rows, axis=0)
