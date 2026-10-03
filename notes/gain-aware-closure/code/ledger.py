# Exact local-defect ledger of the Gaussian closure (telescoping, no linearisation):
#   closure - truth = sum_l (R_l - R_{l+1}),  R_l = closure started from the TRUE moments at layer l,  R_{L+1} = truth.
# Delta_l = R_l - R_{l+1} is the local defect of layer l (Gaussian formula on true moments minus truth),
# carried to the output by the closure itself.
import numpy as np, sys
from math import factorial
from closure import relu_coeffs, relu_var, phi
K = 14
def gstep(mu, S):
    sig = np.sqrt(np.diag(S)); R = S/np.outer(sig, sig); A = relu_coeffs(mu, sig, K)
    C = np.zeros_like(S); Rk = np.ones_like(S)
    for k in range(1, K+1):
        Rk = Rk*R; C += np.outer(A[k], A[k])*Rk/factorial(k)
    np.fill_diagonal(C, relu_var(mu, sig)); return A[0], C
def run_from(Ws, l0, mu, S, keep=False):
    out = []
    for l in range(l0, len(Ws)):
        if l > l0: mu, S = Ws[l] @ m, Ws[l] @ C @ Ws[l].T
        m, C = gstep(mu, S); out.append(m)
    return out if keep else m
if __name__ == "__main__":
    n, L, s = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    Ws = list(np.load(f"W_n{n}_L{L}_s{s}.npy")); tr = np.load(f"truth_n{n}_L{L}_s{s}.npz"); hm = np.load(f"hmom_n{n}_L{L}_s{s}.npz")
    mstar = tr["m"]; Cstar = [hm["S2"][l] - np.outer(hm["S1"][l], hm["S1"][l]) for l in range(L)]
    sstar = [(np.zeros(n), Ws[0] @ Ws[0].T)] + [(Ws[l] @ mstar[l-1], Ws[l] @ Cstar[l-1] @ Ws[l].T) for l in range(1, L)]
    R = [run_from(Ws, l, *sstar[l]) for l in range(L)] + [mstar[-1]]
    mL = mstar[-1]; sc = lambda v: (mL @ v)/(mL @ mL)
    tot = R[0] - mL
    print(f"n={n} L={L} s={s}: closure-truth MSE {np.mean(tot**2):.3e}, scale {sc(tot):+.5f}; telescoping check {np.abs(sum(R[l]-R[l+1] for l in range(L)) - tot).max():.1e}")
    print(" layer | local mean defect: rms  scale(vs m*_l) | carried to output: rms   scale   | split: via mean / via cov (scale)")
    acc = 0
    for l in range(L):
        D = R[l] - R[l+1]; acc += sc(D)
        mloc, Cloc = gstep(*sstar[l]); dm = mloc - mstar[l]
        sl = (mstar[l] @ dm)/(mstar[l] @ mstar[l])
        if l < L-1:
            W = Ws[l+1]
            Dm = run_from(Ws, l+1, W @ mloc, W @ Cstar[l] @ W.T) - R[l+1]
            DC = run_from(Ws, l+1, W @ mstar[l], W @ Cloc @ W.T) - R[l+1]
            split = f"{sc(Dm):+.5f} / {sc(DC):+.5f}"
        else: split = "(final layer)"
        print(f"  {l+1:3d}  |   {np.sqrt(np.mean(dm**2)):.2e}  {sl:+.5f}      |   {np.sqrt(np.mean(D**2)):.2e}  {sc(D):+.5f}  | {split}   cum.scale {acc:+.5f}")
