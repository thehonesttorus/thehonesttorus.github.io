"""Per-neuron variance of point-set probes of the activation means (all layers):

  single     f(x),                                   x ~ N(0, I)
  antithetic (f(x) + f(-x)) / 2                      kills every odd spherical harmonic
  crosspoly  c_n (1/2n) sum_i [f(U e_i) + f(-U e_i)], U Haar   kills degrees 1, 2, 3 exactly (a rotated 3-design)

  python scripts/harmonic_probe.py DATA NET SEED NPAIRS NDESIGNS
Writes $OUT/harm_{NET}_{SEED}.npz with per-estimator sums S1, S2 (L, n) and counts, for pooling across tasks."""
import sys, os, time, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from whest.leaf import c_radial

D, net, seed, npairs, ndes = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
W = [np.ascontiguousarray(w, dtype=np.float32) for w in np.load(f"{D}/W_off{net}.npy")]
L, n = len(W), W[0].shape[0]; rng = np.random.default_rng(10_000 * seed + net); cn = c_radial(n)


def forward(X):                                   # X (n, b) -> activations of every layer, (L, n, b) float64 means only
    out = []
    for Wl in W:
        X = np.maximum(Wl @ X, 0); out.append(X)
    return out


acc = {k: [np.zeros((L, n)), np.zeros((L, n)), 0] for k in ("single", "antithetic", "crosspoly")}
t0 = time.time(); bs = 2048
for _ in range(npairs // bs):
    X = rng.standard_normal((n, bs)).astype(np.float32)
    P, M = forward(X), forward(-X)
    for l in range(L):
        a = P[l].astype(np.float64); b = M[l].astype(np.float64); av = 0.5 * (a + b)
        acc["single"][0][l] += a.sum(1); acc["single"][1][l] += (a * a).sum(1)
        acc["antithetic"][0][l] += av.sum(1); acc["antithetic"][1][l] += (av * av).sum(1)
    acc["single"][2] += bs; acc["antithetic"][2] += bs
t1 = time.time()
for _ in range(ndes):
    U, R = np.linalg.qr(rng.standard_normal((n, n))); U = (U * np.sign(np.diag(R))[None, :]).astype(np.float32)
    P, M = forward(U), forward(-U)
    for l in range(L):
        y = cn * (P[l].astype(np.float64).mean(1) + M[l].astype(np.float64).mean(1)) / 2
        acc["crosspoly"][0][l] += y; acc["crosspoly"][1][l] += y * y
    acc["crosspoly"][2] += 1
print(f"net {net} seed {seed}: {npairs} antithetic pairs in {t1 - t0:.0f}s, {ndes} cross-polytope designs in {time.time() - t1:.0f}s", flush=True)
np.savez(f"{os.environ.get('OUT', '.')}/harm_{net}_{seed}.npz", **{f"{k}_{j}": v[j] for k, v in acc.items() for j in range(3)})
