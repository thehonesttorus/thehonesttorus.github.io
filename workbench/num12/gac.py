# Gain-aware closure (GAC): carry the pre-activations as a scale mixture  z_l = G_l * y_l,
#   y_l ~ N(mu_l, S_l) (conditional Gaussian),  E G_l^2 = 1,  Var G_l^2 = gamma_l,
# with gamma_l = gamma_{l-1} + f_l (critical: ReLU homogeneity carries the gain through unchanged) and
# f_l the weights-only O(n^2) injection of note IX (Prop. 7.4) evaluated on the conditional Gaussian.
# Marginal moments are matched each layer:  E z = E[G] mu,  E zz^T = S + mu mu^T.
# Output mean = E[G_L] * closure mean of y_L.   No fitted constant.
import numpy as np, sys
from math import factorial, lgamma, exp
from closure import relu_coeffs, relu2_coeffs, relu_var, closure
from edgeworth import relu_central
from ledger import gstep
K = 14
def inject(W, mu, S, prev):
    """fresh gain variance created by the source layer prev=(pmu, psig, pR) and seen at pre-activation (mu, S)"""
    n = W.shape[0]; pmu, psig, pR = prev
    q = mu**2 + np.diag(S)
    k3h, k4h = relu_central(pmu, psig)
    W2q = (W*W)/q[:, None]; u = W2q.sum(0); v = (W2q**2).sum(0); t = W.T @ (mu/q)
    A = relu_coeffs(pmu, psig, K); B = relu2_coeffs(pmu, psig, K)
    c = B - 2*A[0]*A; c[0] = 0
    R0 = pR.copy(); np.fill_diagonal(R0, 0.0)
    Cov = np.zeros((n, n)); F = np.zeros((n, n)); Kp = np.zeros((n, n)); Rk = np.ones((n, n))
    for j in range(1, K+1):
        Rk = Rk*R0; f = Rk/factorial(j)
        Cov += np.outer(A[j], A[j])*f; F += np.outer(c[j], c[j])*f; Kp += np.outer(c[j], A[j])*f
    F0 = F - 2*Cov*Cov
    tot = k4h @ (u*u - v) + u @ F0 @ u + 4*np.sum(k3h*u*t) + 4*(u @ Kp @ t)
    return tot/(n*(n-1))
def EG(gam, law="gamma"):
    if gam <= 0: return 1.0
    if law == "first": return 1 - gam/8
    if law == "gamma":          # G^2 ~ Gamma(k, 1/k), Var G^2 = 1/k
        k = 1/gam; return exp(lgamma(k+0.5) - lgamma(k) - 0.5*np.log(k))
    if law == "lognormal":      # log G^2 ~ N(-s/2, s), s = log(1+gam)
        s = np.log1p(gam); return exp(-s/8)
def gac(Ws, gammas=None, law="gamma", tau=1.0):
    """gammas: optional list of gamma_l (pre-activation gain variance per layer) to use instead of the formula"""
    n = Ws[0].shape[0]; out = []; gam = 0.0; g1 = 1.0
    for l, W in enumerate(Ws):
        if l == 0:
            mu = np.zeros(n); S = W @ W.T
        else:
            nu = W @ mbar; Q = W @ Sbar @ W.T                 # marginal mean / second moment of z_l
            if gammas is None:
                mu0 = nu/g1; S0 = Q - np.outer(mu0, mu0)      # conditional state before this layer's injection
                f = inject(W, mu0, S0, prev); gam = tau*gam + f
            else: gam = gammas[l]
            g1 = EG(gam, law)
            mu = nu/g1; S = Q - np.outer(mu, mu)
        M, C = gstep(mu, S)
        mbar = g1*M; Sbar = C + np.outer(M, M)
        sig = np.sqrt(np.diag(S)); prev = (mu, sig, S/np.outer(sig, sig))
        out.append(dict(m=mbar, M=M, gamma=gam, g1=g1, mu=mu, S=S if l == len(Ws)-1 else None))
    return out
if __name__ == "__main__":
    from residue import closure_with_residue
    n, L, s = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    Ws = list(np.load(f"W_n{n}_L{L}_s{s}.npy")); tr = np.load(f"truth_n{n}_L{L}_s{s}.npz"); mt = tr["m"]
    oc = closure(Ws); r1 = closure_with_residue(Ws, tau=1.0); r95 = closure_with_residue(Ws, tau=0.95); g = gac(Ws)
    mse = lambda a, b: np.mean((a-b)**2)
    print(f"n={n} L={L} s={s} (truth T={int(tr['T'])})")
    print(" layer | MSE closure  resid(t=1)  resid(t=.95)  GAC   oracle-scale | 8*scale err: closure  GAC | gamma: formula  GAC")
    for l in range(L):
        mc = oc[l]["m"]; c_or = mc @ (mc-mt[l])/(mc @ mc)
        sc = lambda v: (mt[l] @ (v-mt[l]))/(mt[l] @ mt[l])
        print(f"  {l+1:3d}  | {mse(mc, mt[l]):.3e}  {mse(r1[l]['m'], mt[l]):.3e}  {mse(r95[l]['m'], mt[l]):.3e}  {mse(g[l]['m'], mt[l]):.3e}  {mse(mc*(1-c_or), mt[l]):.3e} |"
              f"  {8*sc(mc):+.4f}  {8*sc(g[l]['m']):+.4f} | {r1[l]['gamma']:.4f}  {g[l]['gamma']:.4f}")
