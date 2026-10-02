"""Monte Carlo truth for the per-layer (2,1) slices D21_l[a,b] = kappa3(z_a, z_a, z_b) of the PRE-activations (row
convention z_l = a_{l-1} @ W_l), in two independent halves (A, B) so that the noise of any comparison can be removed:
||D_est - D_true||^2 is estimated by <D_est - D_A, D_est - D_B> (unbiased, noise-free in expectation).
usage: python3 mc_d21.py MLP LOG2N_PER_HALF OUT.npz"""
import os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '../../fresh-slate/bench'))
import bench
i, lg, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
S = bench.load_set('w1024_d16'); W = bench.weights(S, i).astype(np.float32); L, n, _ = W.shape
Nh = 2 ** lg; ch = 8192
res = {}
t0 = time.time()
for h, seed in (('A', 1001 + i), ('B', 2002 + i)):
    rng = np.random.default_rng(seed)
    s1 = np.zeros((L, n)); s2 = np.zeros((L, n, n)); s21 = np.zeros((L, n, n)); s2d = np.zeros((L, n))
    done = 0
    while done < Nh:
        b = min(ch, Nh - done); a = rng.standard_normal((b, n), dtype=np.float32)
        for l in range(L):
            z = a @ W[l]
            s1[l] += z.sum(0, dtype=np.float64)
            s2[l] += (z.T @ z).astype(np.float64)
            s21[l] += ((z * z).T @ z).astype(np.float64)
            a = np.maximum(z, 0)
        done += b
    m1 = s1 / Nh; M2 = s2 / Nh; M21 = s21 / Nh
    m2d = np.diagonal(M2, axis1=1, axis2=2)
    # kappa3(z_a,z_a,z_b) = E[z_a^2 z_b] - 2 mu_a E[z_a z_b] - mu_b E[z_a^2] + 2 mu_a^2 mu_b
    D = M21 - 2 * m1[:, :, None] * M2 - m2d[:, :, None] * m1[:, None, :] + 2 * (m1 ** 2)[:, :, None] * m1[:, None, :]
    res['D21_' + h] = D.astype(np.float32); res['mu_' + h] = m1; res['var_' + h] = m2d - m1 ** 2
    print('half', h, 'done', round(time.time() - t0, 1), 's', flush=True)
np.savez(out, **res, N_half=Nh)
print('saved', out)
