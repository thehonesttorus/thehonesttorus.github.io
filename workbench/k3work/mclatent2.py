# One Monte Carlo chunk for the transport of the latent fourth cumulant between consecutive layers.
#   python mclatent2.py NET NSAMP CHUNK K LAYERS OUTFILE      (LAYERS = source layers l; each pair is (l, l + 1))
# With mu_l, C_l from mc2_off{NET}_full.npz, U_l the top-K eigenvectors of C_l, Phi_l = Phi(mu_l / sigma_l) and
# W = W_(l+1), per sample:
#   t  = U_l^T d_l                       the latent coordinates at layer l
#   t1 = U_(l+1)^T d_(l+1)               the latent coordinates at layer l + 1 (the truth to be predicted)
#   s1 = U_(l+1)^T W (Phi_l o d_l)       the gated-linear transmission of the whole fluctuation of layer l
#   a1 = U_(l+1)^T W (Phi_l o U_l t)     the transmission of the collective part only, a1 = A^T t
# and the raw moment sums of each K-vector up to fourth order (sum v, sum v v^T, sum v^(x)2 v^T, sum (v^(x)2)(v^(x)2)^T),
# from which latent_an2.py forms the third and fourth cumulant tensors. t1 = s1 + (the nonlinear remainder's read), so
# kappa_4(t1) - kappa_4(s1) is everything born at layer l, and kappa_4(s1) - kappa_4(a1) the idiosyncratic transmission.
import sys, time, numpy as np
net, N, chunk, K = int(sys.argv[1]), int(float(sys.argv[2])), int(sys.argv[3]), int(sys.argv[4])
SRC = [int(x) for x in sys.argv[5].split(",")]; outf = sys.argv[6]
from math import erf, sqrt
Wcol = np.load(f"../official/W_off{net}.npy").astype(np.float32)
L, n, _ = Wcol.shape
T = np.load(f"mc2_off{net}_full.npz")
need = sorted(set(SRC) | {l + 1 for l in SRC})
MU = {l: np.asarray(T["mu"][l], np.float32) for l in need}
U = {}
for l in need:
    ev, V = np.linalg.eigh(np.asarray(T["cov"][l], np.float64)); U[l] = np.ascontiguousarray(V[:, ::-1][:, :K], np.float32)
PHI = {}
for l in SRC:
    s = np.sqrt(np.asarray(T["var"][l], np.float64)); a = np.asarray(T["mu"][l], np.float64) / s
    PHI[l] = np.array([0.5 * (1 + erf(x / sqrt(2))) for x in a], np.float32)
G = {l: np.ascontiguousarray(U[l + 1].T @ (Wcol[l + 1] * PHI[l][None, :]), np.float32) for l in SRC}   # K x n
A = {l: np.ascontiguousarray(G[l] @ U[l], np.float32) for l in SRC}                                     # K x K
iu = np.triu_indices(K)
names = ("t", "t1", "s1", "a1")
acc = {(l, v): dict(m1=np.zeros(K), m2=np.zeros((K, K)), m3=np.zeros((len(iu[0]), K)), m4=np.zeros((len(iu[0]), len(iu[0]))))
       for l in SRC for v in names}
rng = np.random.default_rng([4712, net, chunk])
B = 16384; done = 0; t0 = time.time()


def add(a, X):
    # X: K x b. Pairs p <= q only (the symmetric square), so the fourth moments are a (K(K+1)/2)^2 matrix.
    P = (X[iu[0]] * X[iu[1]])                      # K(K+1)/2 x b
    a["m1"] += X.sum(1, dtype=np.float64); a["m2"] += X @ X.T
    a["m3"] += P @ X.T; a["m4"] += P @ P.T


while done < N:
    b = min(B, N - done)
    Y = rng.standard_normal((n, b), dtype=np.float32)
    D = {}
    for l in range(max(need) + 1):
        Z = Wcol[l] @ Y
        if l in MU:
            D[l] = Z - MU[l][:, None]
        Y = np.maximum(Z, 0.0, out=Z)
    for l in SRC:
        t = U[l].T @ D[l]; t1 = U[l + 1].T @ D[l + 1]
        s1 = G[l] @ D[l]; a1 = A[l] @ t
        for v, X in (("t", t), ("t1", t1), ("s1", s1), ("a1", a1)):
            add(acc[(l, v)], X)
    done += b
    if (done // B) % 8 == 0:
        print(f"{done} samples, {time.time() - t0:.0f}s", flush=True)
np.savez(outf, c=done, src=np.array(SRC), K=K, **{f"A_{l}": A[l] for l in SRC},
         **{f"{k}_{v}_{l}": acc[(l, v)][k] for l in SRC for v in names for k in acc[(l, v)]})
print(f"done {done} samples in {time.time() - t0:.0f}s", flush=True)
