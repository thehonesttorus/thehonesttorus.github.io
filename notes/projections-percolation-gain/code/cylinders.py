# Test of the "cylinder mixture" (tiling) proposal: condition on sign(Px) for k frame directions,
# run the Gaussian closure inside each cylinder (truncated-Gaussian moments), average.
import numpy as np, sys, time
from math import factorial
from closure import relu_coeffs, relu_var
def closure_from(Ws, m0, S0, K=14):
    n = Ws[0].shape[0]; mu = Ws[0] @ m0; S = Ws[0] @ S0 @ Ws[0].T
    for l, W in enumerate(Ws):
        if l > 0: mu = W @ m; S = W @ C @ W.T
        sig = np.sqrt(np.diag(S)); R = S/np.outer(sig, sig)
        A = relu_coeffs(mu, sig, K); m = A[0]
        C = np.zeros((n, n)); Rk = np.ones((n, n))
        for k in range(1, K+1):
            Rk = Rk*R; C += np.outer(A[k], A[k])*Rk/factorial(k)
        np.fill_diagonal(C, relu_var(mu, sig))
    return m
def frame(Ws, kind, k):
    n = Ws[0].shape[0]
    if kind == "W1":      # top right singular vectors of W1
        _, _, Vt = np.linalg.svd(Ws[0]); return Vt[:k]
    if kind == "jac":     # top right singular vectors of the mean Jacobian (1/2 W_L)...(1/2 W_1)
        J = np.eye(n)
        for W in Ws: J = 0.5*W @ J
        _, _, Vt = np.linalg.svd(J); return Vt[:k]
    if kind == "rand":
        Q, _ = np.linalg.qr(np.random.default_rng(7).standard_normal((n, k))); return Q.T
def mixture(Ws, P):
    k, n = P.shape; c = np.sqrt(2/np.pi)
    S0 = np.eye(n) - (2/np.pi)*P.T @ P
    out = 0
    for idx in range(2**k):
        eps = np.array([1.0 if (idx >> j) & 1 else -1.0 for j in range(k)])
        out = out + closure_from(Ws, P.T @ (c*eps), S0)/2**k
    return out
if __name__ == "__main__":
    n, L, s = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    Ws = list(np.load(f"W_n{n}_L{L}_s{s}.npy")); tr = np.load(f"truth_n{n}_L{L}_s{s}.npz"); m = tr["m"][-1]
    base = closure_from(Ws, np.zeros(n), np.eye(n))
    print(f"n={n} L={L} s={s}: MC@B {tr['v'][-1].mean()/65536:.2e}  closure {np.mean((base-m)**2):.3e}", flush=True)
    for kind in ["W1", "jac", "rand"]:
        for k in [2, 4, 6]:
            t0 = time.time(); mm = mixture(Ws, frame(Ws, kind, k))
            print(f"   frame {kind:4s} k={k}: mixture MSE {np.mean((mm-m)**2):.3e}   ({time.time()-t0:.0f}s)", flush=True)
