# One chunk of the Monte Carlo of mcstats.py, for parallel runs (merged by mcmerge.py into the files mcstats.py writes).
#   python mcchunk.py NET NSAMP SEED CHUNK OUTFILE [post]
# Samples come from default_rng([SEED, CHUNK]), so chunks are independent. Every chunk uses the same shift: the
# dataset's own all-layer means m_l for y_l and W_l m_(l-1) for z_l, computed in float64 and rounded to multiples of
# 2^-10. Rounding makes the shift bitwise identical on every host whatever its BLAS, so the raw-moment sums of
# d = z - m0 add exactly across chunks. The accumulators are those of mcstats.py: s1..s4, m11, m21, m22, m31, yfin
# (and the _y set with post).
import sys, time, numpy as np
net, N, seed, chunk, outf = int(sys.argv[1]), int(float(sys.argv[2])), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
POST = len(sys.argv) > 6 and sys.argv[6] == "post"
Wcol = np.load(f"../official/W_off{net}.npy").astype(np.float32)
L, n, _ = Wcol.shape
B = 32768
KEYS = ("s1", "s2", "s3", "s4", "m11", "m21", "m22", "m31")


def newacc():
    return dict(c=0, s1=np.zeros((L, n)), s2=np.zeros((L, n)), s3=np.zeros((L, n)), s4=np.zeros((L, n)),
                m11=np.zeros((L, n, n)), m21=np.zeros((L, n, n)), m22=np.zeros((L, n, n)), m31=np.zeros((L, n, n)), yfin=np.zeros(n))


mt = np.load(f"../official/truth_off{net}.npz")["m"].astype(np.float64)
grid = lambda v: (np.round(v * 1024.0) / 1024.0).astype(np.float32)
m0y = grid(mt)
m0 = grid(np.stack([np.zeros(n)] + [Wcol[l].astype(np.float64) @ mt[l - 1] for l in range(1, L)]))
rng = np.random.default_rng([seed, chunk])
a, ay = newacc(), (newacc() if POST else None)
t0 = time.time(); done = 0; nb = 0
while done < N:
    b = min(B, N - done)
    Y = rng.standard_normal((n, b), dtype=np.float32)
    a["c"] += b
    if POST:
        ay["c"] += b
    for l in range(L):
        Z = Wcol[l] @ Y
        D = Z - m0[l][:, None]; D2 = D * D
        a["s1"][l] += D.sum(1, dtype=np.float64); a["s2"][l] += D2.sum(1, dtype=np.float64)
        a["s3"][l] += (D2 * D).sum(1, dtype=np.float64); a["s4"][l] += (D2 * D2).sum(1, dtype=np.float64)
        a["m11"][l] += D @ D.T; a["m21"][l] += D2 @ D.T; a["m22"][l] += D2 @ D2.T; a["m31"][l] += (D2 * D) @ D.T
        Y = np.maximum(Z, 0.0, out=Z)
        if POST:
            Dy = Y - m0y[l][:, None]; Dy2 = Dy * Dy
            ay["s1"][l] += Dy.sum(1, dtype=np.float64); ay["s2"][l] += Dy2.sum(1, dtype=np.float64)
            ay["s3"][l] += (Dy2 * Dy).sum(1, dtype=np.float64); ay["s4"][l] += (Dy2 * Dy2).sum(1, dtype=np.float64)
            ay["m11"][l] += Dy @ Dy.T; ay["m21"][l] += Dy2 @ Dy.T; ay["m22"][l] += Dy2 @ Dy2.T; ay["m31"][l] += (Dy2 * Dy) @ Dy.T
    a["yfin"] += Y.sum(1, dtype=np.float64)
    done += b; nb += 1
    if nb % 16 == 0:
        print(f"{done} samples, {time.time() - t0:.0f}s", flush=True)
save = dict(c=a["c"], m0=m0, m0y=m0y, yfin=a["yfin"], **{k: a[k] for k in KEYS})
if POST:
    save.update(cy=ay["c"], **{f"{k}_y": ay[k] for k in KEYS})
np.savez(outf, **save)
print(f"chunk {chunk}: {N} samples in {time.time() - t0:.0f}s", flush=True)
