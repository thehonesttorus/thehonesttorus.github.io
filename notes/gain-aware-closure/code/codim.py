# Codimension test of the "dual truncation" proposal.  One closure step z = W g(y), y Gaussian (closure state):
#   single-wall terms (codim 1): D3, D4 (one source neuron's kink)      pair-of-walls (codim 2): P3, T3, PP4, PD4
# (a) Edgeworth mean shift from each group: rms and scale part (8 x projection on the closure mean)
# (b) the coherent gain injection f_l split the same way: k4-diagonal and k3mu (one wall) vs F0 and Kp (two walls)
# (c) "only O(n) of the O(n^2) pairs matter": share of the pair sum u^T F0 u carried by the top n / 10n pairs,
#     ranked by |contribution| and by |correlation|
import numpy as np, sys
from math import factorial
from closure import closure, relu_coeffs, relu2_coeffs
from edgeworth import cumulants_next, edgeworth_shift, relu_central
n, L, s = int(sys.argv[1]), 16, int(sys.argv[2])
Ws = list(np.load(f"W_n{n}_L{L}_s{s}.npy")); oc = closure(Ws, keep=True); K = 14
print(f"n={n} s={s}")
print(" layer | shift rms: codim1  codim2 | 8*scale: codim1  codim2 | gain f_l: 1-wall  2-wall | pair sum share: top n  top 10n (by contribution) | top n  top 10n (by |rho|)")
for l in [1, 3, 7, 11, 15]:
    p = oc[l-1]; W = Ws[l]; mu, sig, m = oc[l]["mu"], oc[l]["sig"], oc[l]["m"]
    _, _, out = cumulants_next(W, p["mu"], p["sig"], p["R"], J=2)
    s1 = -edgeworth_shift(mu, sig, out["D3"], out["D4"]); s2 = -edgeworth_shift(mu, sig, out["P3"]+out["T3"], out["PP4"]+out["PD4"])
    sc = lambda v: 8*(m @ v)/(m @ m)
    # gain injection split (as in residue.py / gac.inject)
    pmu, psig, pR = p["mu"], p["sig"], p["R"]; q = mu**2 + sig**2
    k3h, k4h = relu_central(pmu, psig)
    W2q = (W*W)/q[:, None]; u = W2q.sum(0); v = (W2q**2).sum(0); t = W.T @ (mu/q)
    A = relu_coeffs(pmu, psig, K); B = relu2_coeffs(pmu, psig, K); c = B - 2*A[0]*A; c[0] = 0
    R0 = pR.copy(); np.fill_diagonal(R0, 0.0)
    Cov = np.zeros((n, n)); F = np.zeros((n, n)); Kp = np.zeros((n, n)); Rk = np.ones((n, n))
    for j in range(1, K+1):
        Rk = Rk*R0; f = Rk/factorial(j); Cov += np.outer(A[j], A[j])*f; F += np.outer(c[j], c[j])*f; Kp += np.outer(c[j], A[j])*f
    F0 = F - 2*Cov*Cov
    one = k4h @ (u*u - v) + 4*np.sum(k3h*u*t); X = np.outer(u, u)*F0; two = X.sum() + 4*(u @ Kp @ t)
    iu = np.triu_indices(n, 1); x = 2*X[iu]; tot = x.sum(); r = np.abs(R0[iu])
    o1 = np.argsort(-np.abs(x)); o2 = np.argsort(-r)
    share = lambda o, k: x[o[:k]].sum()/tot
    print(f"  {l+1:3d}  |   {np.sqrt(np.mean(s1**2)):.1e}  {np.sqrt(np.mean(s2**2)):.1e} |  {sc(s1):+.4f}  {sc(s2):+.4f} |  {one/(n*n-n):+.5f}  {two/(n*n-n):+.5f} |"
          f"  {share(o1, n):.3f}  {share(o1, 10*n):.3f}                     |  {share(o2, n):.3f}  {share(o2, 10*n):.3f}", flush=True)
