# GAC + fresh kappa_3 trees at reduced cost: D3 (O(n^2)), P3 with Hermite orders j <= jmax (jmax matmuls),
# T3 at J = 1 (2 matmuls).  Mean shift only (Edgeworth, kappa_3 and kappa_3^2 terms).
import numpy as np, sys
from math import factorial
from ledger import gstep
from gac import inject, EG
from closure import relu_coeffs, relu2_coeffs
from edgeworth import relu_central, edgeworth_shift
def k3_cheap(W, mu, sig, R, jmax=2, parts=("D3", "P3", "T3")):
    A = relu_coeffs(mu, sig, max(jmax, 2)); B = relu2_coeffs(mu, sig, max(jmax, 2)); c = B - 2*A[0]*A
    R0 = R.copy(); np.fill_diagonal(R0, 0.0); W2 = W*W; k3h, _ = relu_central(mu, sig); k3 = np.zeros(len(mu))
    if "D3" in parts: k3 += (W2*W) @ k3h
    if "P3" in parts:
        Rj = np.ones_like(R0)
        for j in range(1, jmax+1):
            Rj = Rj*R0; k3 += 3*np.sum(((W2*c[j]) @ (Rj/factorial(j)))*(W*A[j]), 1)
    if "T3" in parts:
        U1 = (W*A[1]) @ R0; V11 = (W2*A[1]**2) @ (R0*R0); k3 += 3*np.sum(W*A[2]*(U1*U1 - V11), 1)
    return k3
def gac_k3(Ws, jmax=2, parts=("D3", "P3", "T3"), lmax=99):
    n = Ws[0].shape[0]; out = []; gam = 0.0; g1 = 1.0
    for l, W in enumerate(Ws):
        if l == 0: mu = np.zeros(n); S = W @ W.T
        else:
            nu = W @ mbar; Q = W @ Sbar @ W.T
            mu0 = nu/g1; S0 = Q - np.outer(mu0, mu0); gam = gam + inject(W, mu0, S0, prev)
            g1 = EG(gam); mu = nu/g1; S = Q - np.outer(mu, mu)
        M, C = gstep(mu, S); sig = np.sqrt(np.diag(S))
        if 0 < l <= lmax:
            k3 = k3_cheap(W, *prev, jmax=jmax, parts=parts); M = M + edgeworth_shift(mu, sig, k3, np.zeros(n))
        mbar = g1*M; Sbar = C + np.outer(M, M)
        prev = (mu, sig, S/np.outer(sig, sig))
        out.append(dict(m=mbar, gamma=gam))
    return out
if __name__ == "__main__":
    for tag in sys.argv[1].split(","):
        Ws = list(np.load(f"W_{tag}.npy").astype(np.float64)); mt = np.load(f"truth_{tag}.npz")["m"][-1].astype(np.float64)
        mse = lambda a: np.mean((a-mt)**2); row = f"{tag}:"
        for kw in [dict(jmax=1), dict(jmax=2), dict(jmax=4), dict(jmax=2, lmax=6), dict(jmax=1, parts=("D3", "P3"))]:
            row += f" | {kw} {mse(gac_k3(Ws, **kw)[-1]['m']):.3e}"
        print(row, flush=True)
