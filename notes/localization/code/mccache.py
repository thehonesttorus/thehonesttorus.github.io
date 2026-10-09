# Monte Carlo statistics cache for one network (offline; plain numpy, float32 forward, float64 accumulation).
#   python mccache.py NET [LOG2_SAMPLES=20] [BATCH=4096]
# Per layer l = 0..15, with h_l the pre-activation and x_{l+1} = relu(h_l) (truth 'm'[l] is E x_{l+1}):
#   s1..s4  : per-neuron raw power sums of h_l (float64)        -> mean, var, kappa3, kappa4
#   pos     : per-neuron count of h_l > 0                       -> gate probability
#   xs1     : per-neuron sum of x_{l+1}                          -> check against truth
#   G       : E[X^T x_{l+1}] (n x n, input index first)          -> first-chaos (Stein) map E[grad x_{l+1}]
#   Hc      : E[h_l^T h_l] (n x n)                               -> pre-activation second moments / covariance
#   rad     : E[(|X|^2 - n) x_{l+1}] / sqrt(2n)                  -> radial (trace) part of the second chaos
# Saved to $OUT/mccache_{NET}.npz (n x n blocks as float32).
import os, sys, time
import numpy as np

net = int(sys.argv[1]); lg = int(sys.argv[2]) if len(sys.argv) > 2 else 20
bs = int(sys.argv[3]) if len(sys.argv) > 3 else 4096
S = 1 << lg; nb = S // bs
Wcol = np.load(f"../official/W_off{net}.npy")
ws = [np.ascontiguousarray(W.T).astype(np.float32) for W in Wcol]
L, n = len(ws), ws[0].shape[0]
s = {k: np.zeros((L, n)) for k in ("s1", "s2", "s3", "s4", "pos", "xs1", "rad")}
G = np.zeros((L, n, n)); Hc = np.zeros((L, n, n))
rng = np.random.default_rng(1000003 * (net + 1))
t0 = time.time()
for b in range(nb):
    X = rng.standard_normal((bs, n), dtype=np.float32)
    r = ((X.astype(np.float64) ** 2).sum(1) - n) / np.sqrt(2.0 * n)
    x = X
    for l in range(L):
        h = x @ ws[l]
        hd = h.astype(np.float64); h2 = hd * hd
        s["s1"][l] += hd.sum(0); s["s2"][l] += h2.sum(0); s["s3"][l] += (h2 * hd).sum(0); s["s4"][l] += (h2 * h2).sum(0)
        s["pos"][l] += (h > 0).sum(0)
        x = np.maximum(h, np.float32(0.0))
        s["xs1"][l] += x.sum(0, dtype=np.float64)
        s["rad"][l] += r @ x.astype(np.float64)
        G[l] += (X.T @ x).astype(np.float64)
        Hc[l] += (h.T @ h).astype(np.float64)
    if b in (0, nb // 4, nb // 2, nb - 1):
        print(f"net {net}: batch {b + 1}/{nb} at {time.time() - t0:.0f}s", flush=True)
mt = np.load(f"../official/truth_off{net}.npz")["m"].astype(np.float64)
xm = s["xs1"] / S
print(f"net {net}: S=2^{lg}, MSE of MC means vs truth per layer (expect var/S): "
      + " ".join(f"{np.mean((xm[l] - mt[l]) ** 2):.1e}" for l in (0, 5, 10, 15)), flush=True)
out = os.environ.get("OUT", ".")
np.savez(f"{out}/mccache_{net}.npz", S=S, **s, G=(G / S).astype(np.float32), Hc=(Hc / S).astype(np.float32))
print(f"net {net}: saved in {time.time() - t0:.0f}s", flush=True)
