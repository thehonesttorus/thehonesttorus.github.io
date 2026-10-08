# The proposal's fibre is a rank-32 covariance. What does truncating the closure covariance do at width 1024?
#   rank-r      : C <- top-r eigen-part of C (the proposal)
#   diag+rank-r : C <- top-r eigen-part + exact diagonal (generous)
import numpy as np, sys
from ledger import gstep
n, L, s = int(sys.argv[1]), 16, int(sys.argv[2])
Ws = list(np.load(f"W_n{n}_L{L}_s{s}.npy")); mt = np.load(f"truth_n{n}_L{L}_s{s}.npz")["m"][-1]
def run(r=None, keepdiag=False):
    mu = np.zeros(n); S = Ws[0] @ Ws[0].T
    for l in range(L):
        if l > 0: mu, S = Ws[l] @ m, Ws[l] @ C @ Ws[l].T
        m, C = gstep(mu, S)
        if r is not None:
            ev, V = np.linalg.eigh(C); Cr = (V[:, -r:]*ev[-r:]) @ V[:, -r:].T
            if keepdiag: Cr[np.diag_indices(n)] = np.diag(C)
            C = Cr
    return m
mse = lambda a: np.mean((a-mt)**2)
print(f"n={n} s={s}: full closure MSE {mse(run()):.3e}")
for r in [32, 128]:
    print(f"  rank-{r}: MSE {mse(run(r)):.3e}    diag+rank-{r}: MSE {mse(run(r, True)):.3e}", flush=True)
