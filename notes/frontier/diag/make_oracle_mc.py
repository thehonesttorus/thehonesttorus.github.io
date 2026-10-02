"""Plumbing/low-precision oracle: per-layer marginals (mu, var, k3, k4) of the PRE-activations by plain Monte Carlo
(row-vector convention z_l = a_{l-1} @ W_l).  N must be huge for real diagnostics (use the N=1e9 HF moments instead).
usage: python3 make_oracle_mc.py MLP N OUT.npz"""
import sys, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '../../fresh-slate/bench'))
import bench
i, N, out = int(sys.argv[1]), int(float(sys.argv[2])), sys.argv[3]
S = bench.load_set('w1024_d16'); W = bench.weights(S, i).astype(np.float32); L, n, _ = W.shape
s1 = np.zeros((L, n)); s2 = np.zeros((L, n)); s3 = np.zeros((L, n)); s4 = np.zeros((L, n))
rng = np.random.default_rng(12345 + i); ch = 8192; done = 0
while done < N:
    b = min(ch, N - done); a = rng.standard_normal((b, n), dtype=np.float32)
    for l in range(L):
        z = a @ W[l]; zd = z.astype(np.float64)
        s1[l] += zd.sum(0); z2 = zd * zd; s2[l] += z2.sum(0); s3[l] += (z2 * zd).sum(0); s4[l] += (z2 * z2).sum(0)
        a = np.maximum(z, 0)
    done += b
m1, m2, m3, m4 = s1 / N, s2 / N, s3 / N, s4 / N
var = m2 - m1 ** 2
k3 = m3 - 3 * m1 * m2 + 2 * m1 ** 3
k4 = m4 - 4 * m1 * m3 - 3 * m2 ** 2 + 12 * m1 ** 2 * m2 - 6 * m1 ** 4
np.savez(out, mu=m1, var=var, k3=k3, k4=k4, N=N)
print('saved', out, 'N', N)
