# Closure round, sub-ledger from the v29 dump (layers 5, 10, 14, network 0):
#  (a) the chain's fitted fourth-cumulant sector (g4row = k4 diag, wk4m = K22 slice, wk431 = K31 slice of the pre-activation)
#      against the derived scale-mixture sector of E10 (z = G y, Var G^2 = g, g read from the chain's own D3), all slices;
#  (b) the size of the Mehler-3 and Mehler-4 terms the pair programs drop, for the (1,1) and (2,1) kernels.
import numpy as np, pickle, sys, math
sys.path.insert(0, "../ncgprob"); sys.path.insert(0, "../num12")
from pairvar import sm_slices
from closure import relu_coeffs, relu2_coeffs
from scipy.special import ndtr
phi = lambda x: np.exp(-0.5 * x * x) / np.sqrt(2 * np.pi)
net = int(sys.argv[1]) if len(sys.argv) > 1 else 0
D = pickle.load(open(f"v29dump_off{net}.pkl", "rb"))
cos = lambda A, B: float(np.sum(A * B) / np.sqrt(np.sum(A * A) * np.sum(B * B)))
nrm = np.linalg.norm
for d in D:
    l = d["layer"]; mu, var, C_off, D3 = d["mu"], d["var"], d["C_off"], d["D3"]; n = len(mu)
    sig = np.sqrt(var); a = mu / sig; Ph = ndtr(a); ph = phi(a); m = sig * (a * Ph + ph)
    # g as the chain's hook defines it: D3 ~ 1.5 g mu var
    v = 1.5 * mu * var; g = float(v @ D3 / (v @ v)); corr = cos(v, D3)
    S = C_off + np.diag(var)
    sm = sm_slices(mu, S, g)
    k4c, K22c, K31c = d["g4row"], d["wk4m"], d.get("wk431", None)
    k4d, K22d, K31d = sm["k4"], sm["K22"].copy(), sm["K31"].copy(); np.fill_diagonal(K22d, 0); np.fill_diagonal(K31d, 0)
    print(f"\nlayer {l}: g = {g:.4f} (corr D3 vs 1.5 mu var {corr:.3f}); |mu| {nrm(mu):.3f} mean alpha {a.mean():.3f}")
    print(f"  k4 diag: chain/derived norm ratio {nrm(k4c)/nrm(k4d):.3f}, cosine {cos(k4c, k4d):.4f}; chain/3 g sig^4 ratio {nrm(k4c)/nrm(3*g*var*var):.3f}")
    print(f"  K22 off: chain/derived norm ratio {nrm(K22c)/nrm(K22d):.3f}, cosine {cos(K22c, K22d):.4f}; derived vs g sig^2 sig^2^T cosine {cos(K22d, np.outer(var, var) - np.diag(var*var)):.4f}, vs C_off*C_off {cos(K22d, C_off*C_off):.4f}")
    # mean-coupled part of the derived K22: regress on the rank-one sig2 sig2^T and on mu-structures
    if K31c is not None: print(f"  K31 off: chain/derived norm ratio {nrm(K31c)/nrm(K31d):.3f}, cosine {cos(K31c, K31d):.4f}")
    if True:
        print(f"  K31 off derived: |K31d|/|K22d| {nrm(K31d)/nrm(K22d):.3f}; vs C_off cosine {cos(K31d, C_off):.4f}, vs d(var) C_off {cos(K31d, var[:, None]*C_off):.4f}, vs d(mu) C_off d(mu) {cos(K31d, (mu[:, None]*C_off)*mu[None, :]):.4f}")
    # (b) Mehler terms of the (1,1) and (2,1) kernels in powers of C_off (coefficients: shifted Hermite of relu, relu^2)
    A = relu_coeffs(mu, sig, 5); B = relu2_coeffs(mu, sig, 5)   # A[k] = E[relu He_k], includes sig^k? -> coefficient of rho^k is A_k A_k / k! with rho = C/(sig sig)
    R = C_off / np.outer(sig, sig)
    def term(P, Q, k): return np.outer(P[k], Q[k]) * (R ** k) / math.factorial(k)
    T11 = [term(A, A, k) for k in (1, 2, 3, 4)]; T21 = [term(B, A, k) for k in (1, 2, 3, 4)]
    rho = np.abs(R[~np.eye(n, dtype=bool)]); print(f"  |rho| off-diagonal: rms {np.sqrt(np.mean(rho**2)):.4f}, 99.9% {np.quantile(rho, 0.999):.4f}, max {rho.max():.4f}")
    for name, T in (("(1,1)", T11), ("(2,1)", T21)):
        n1 = nrm(T[0]); print(f"  {name} Mehler terms |k|/|k=1|: k=2 {nrm(T[1])/n1:.2e}, k=3 {nrm(T[2])/n1:.2e}, k=4 {nrm(T[3])/n1:.2e}; coherent (sum over entries) k=3/k=1 {T[2].sum()/T[0].sum():+.2e}, k=4/k=1 {T[3].sum()/T[0].sum():+.2e}")
