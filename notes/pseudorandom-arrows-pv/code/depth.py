# Depth: reduce layer l only (keep sketch W P with P = proj onto mean direction, plus row norms),
# replace the residual by (a) fresh Gaussian, (b) one shared direction, (c) rows in a random d-dim subspace.
# Compare the observed final-layer MSE with the annealed-downstream identity (NNGP kernel form).
import numpy as np, sys
spec = np.load("spec.npy")                    # spec[D-1] = Taylor coeffs of rho^{oD}
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
n, L, l, N = 256, 8, 3, 60000
s2 = 2.0 / n
Ws = [rng.standard_normal((n, n)) * np.sqrt(s2) for _ in range(L)]
X = rng.standard_normal((N, n))
def fwd(H, layers):
    for W in layers: H = np.maximum(H @ W.T, 0)
    return H
Hlm1 = fwd(X, Ws[:l-1])                       # h_{l-1}
m = Hlm1.mean(0); e = m / np.linalg.norm(m); P = np.outer(e, e); Q = np.eye(n) - P
W = Ws[l-1]; rho = np.linalg.norm(W @ Q, axis=1)
def rowsph(Ures):                              # rescale residual rows to the true residual norms rho_j
    return Ures * (rho / np.linalg.norm(Ures, axis=1))[:, None]
red = {}
red["fresh Gaussian residual"] = W @ P + rowsph(rng.standard_normal((n, n)) @ Q)
g = Q @ rng.standard_normal(n)
red["shared direction (rank 2)"] = W @ P + rowsph(np.tile(g, (n, 1)))
for d in [4, 32]:
    B = Q @ rng.standard_normal((n, d))
    red[f"rows in random {d}-dim subspace"] = W @ P + rowsph(rng.standard_normal((n, d)) @ B.T)
D = L - l
c = np.polynomial.polynomial
def kern(A, B):                               # |a||b| rho^{oD}(cos)  (NNGP limit of the downstream two-point function)
    na = np.linalg.norm(A, axis=1); nb = np.linalg.norm(B, axis=1)
    C = np.clip((A @ B.T) / np.maximum(np.outer(na, nb), 1e-300), -1, 1)
    return np.outer(na, nb) * c.polyval(C, spec[D-1])
Htrue = np.maximum(Hlm1 @ W.T, 0)
fT = fwd(Htrue, Ws[l:]).mean(0)
idx = rng.choice(N, 1500, replace=False)
KTT = kern(Htrue[idx], Htrue[idx]); np.fill_diagonal(KTT, 0)
print(f"n={n}, L={L}, reduced layer l={l} (downstream D={D});  var of f_j across neurons = {fT.var():.3e}")
for name, Wr in red.items():
    Hr = np.maximum(Hlm1 @ Wr.T, 0)
    fR = fwd(Hr, Ws[l:]).mean(0)
    mse = np.mean((fT - fR)**2)
    KRR = kern(Hr[idx], Hr[idx]); np.fill_diagonal(KRR, 0)
    KTR = kern(Htrue[idx], Hr[idx[::-1]])     # independent pairs (X, X')
    k = len(idx)
    pred = (KTT.sum() / (k*(k-1)) + KRR.sum() / (k*(k-1)) - 2 * KTR.mean()) / n
    print(f"  {name:32s} observed MSE = {mse:.3e}   identity (NNGP kernel) prediction = {pred:.3e}")
