# The proposed width-64, depth-8 probe of the "one filtration" picture.
#   closure        : Gaussian closure (codimension-0 page)
#   wall page      : closure + the diagonal wall current (single-neuron kappa_3, kappa_4 one-point corrections, D3 + D4) at
#                    every layer -- the conjecture's tau^(1)
#   cylinder page  : 16 cylinders = sign patterns of 4 layer-1 neurons; exact cell-conditional moments of h_1 (4e6 samples),
#                    closure from layer 2 in each cell, mixed with the cell weights
#   GAC            : gain-aware closure (note XI)
#   d1             : wall page - cylinder page
# Question: does the closure residual (truth - closure) correlate with d1, and which page captures it?
import numpy as np, sys
from ledger import gstep
from edgeworth import cumulants_next, edgeworth_shift
from gac import gac
from scipy.special import ndtr
n, L = 64, 8
def pages(Ws, mu, S, wall=False, l0=0):
    for l in range(l0, len(Ws)):
        if l > l0: mu, S = Ws[l] @ m, Ws[l] @ C @ Ws[l].T
        sig = np.sqrt(np.diag(S)); R = S/np.outer(sig, sig)
        m, C = gstep(mu, S)
        if wall and l > l0:
            _, _, out = cumulants_next(Ws[l], pmu, psig, pR, J=2, terms=("D3", "D4"))
            m = m + edgeworth_shift(mu, sig, out["D3"], out["D4"])
        pmu, psig, pR = mu, sig, R
    return m
corr = lambda a, b: np.corrcoef(a, b)[0, 1]
print(" seed | MSE: closure   wall      cylinder  GAC      | corr with residual: d1    wall-cl  cyl-cl  GAC-cl | min_l 2<P(1-P)>")
for s in [int(x) for x in sys.argv[1].split(",")]:
    Ws = list(np.load(f"W_n{n}_L{L}_s{s}.npy")); mt = np.load(f"truth_n{n}_L{L}_s{s}.npz")["m"][-1]
    mc = pages(Ws, np.zeros(n), Ws[0] @ Ws[0].T); mw = pages(Ws, np.zeros(n), Ws[0] @ Ws[0].T, wall=True)
    rng = np.random.default_rng(s + 7); mcyl = np.zeros(n); N = 0; cells = {}
    for _ in range(40):
        x = rng.standard_normal((100000, n)); h = np.maximum(x @ Ws[0].T, 0); lab = (h[:, :4] > 0) @ (2**np.arange(4))
        for c in range(16):
            X = h[lab == c]; a = cells.setdefault(c, [0, np.zeros(n), np.zeros((n, n))]); a[0] += len(X); a[1] += X.sum(0); a[2] += X.T @ X
    N = sum(a[0] for a in cells.values())
    for c, (k, s1, s2) in cells.items():
        m1 = s1/k; C1 = s2/k - np.outer(m1, m1)
        mcyl += (k/N)*pages(Ws, Ws[1] @ m1, Ws[1] @ C1 @ Ws[1].T, l0=1)
    mg = gac(Ws)[-1]["m"]
    # non-scale forgetting rate per layer (note IX): 2 <P(1-P)> with P = Phi(mu/sigma) along the closure
    mu, S = np.zeros(n), Ws[0] @ Ws[0].T; rates = []
    for l in range(L):
        if l > 0: mu, S = Ws[l] @ m_, Ws[l] @ C_ @ Ws[l].T
        P = ndtr(mu/np.sqrt(np.diag(S))); rates.append(2*np.mean(P*(1-P))); m_, C_ = gstep(mu, S)
    res = mt - mc; mse = lambda a: np.mean((a-mt)**2)
    print(f" {s:4d} | {mse(mc):.2e}  {mse(mw):.2e}  {mse(mcyl):.2e}  {mse(mg):.2e} |  {corr(mw-mcyl, res):+.2f}   {corr(mw-mc, res):+.2f}   {corr(mcyl-mc, res):+.2f}   {corr(mg-mc, res):+.2f}  | {min(rates[1:]):.3f}", flush=True)
