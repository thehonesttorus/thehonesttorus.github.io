# Impulse response of the closure: how a perturbation of the post-activation means at layer l reaches the
# output scale. Uniform (m_l -> (1+eps) m_l) and gap-shaped (delta m = eps*sigma*phi(a)(1+a^2)) perturbations.
import numpy as np, sys
from ledger import gstep, run_from
from closure import phi
n, L, s = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
Ws = list(np.load(f"W_n{n}_L{L}_s{s}.npy"))
# closure trajectory
mus, Ss, ms, Cs = [], [], [], []
mu, S = np.zeros(n), Ws[0] @ Ws[0].T
for l in range(L):
    if l > 0: mu, S = Ws[l] @ m, Ws[l] @ C @ Ws[l].T
    m, C = gstep(mu, S); mus.append(mu); Ss.append(S); ms.append(m); Cs.append(C)
mL = ms[-1]; sc = lambda v: (mL @ v)/(mL @ mL)
eps = 1e-3
print("layer   c_l    t_uniform   t_gap(per unit scale-content)   [output scale per unit local scale]")
for l in range(L-1):
    W = Ws[l+1]; sig = np.sqrt(np.diag(Ss[l])); a = mus[l]/sig
    cfrac = (ms[l] @ ms[l])/np.sum(ms[l]**2 + np.diag(Cs[l]))
    base = run_from(Ws, l+1, W @ ms[l], W @ Cs[l] @ W.T)
    up = run_from(Ws, l+1, W @ (ms[l]*(1+eps)), W @ Cs[l] @ W.T)
    g = sig*phi(a)*(1+a*a); g_sc = (ms[l] @ g)/(ms[l] @ ms[l])
    ug = run_from(Ws, l+1, W @ (ms[l] + eps*g/g_sc), W @ Cs[l] @ W.T)
    print(f"{l+1:4d}  {cfrac:.3f}   {sc(up-base)/eps:+.3f}      {sc(ug-base)/eps:+.3f}")
