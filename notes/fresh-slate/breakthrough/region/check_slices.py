import sys, numpy as np
sys.path.insert(0, '../../bench')
import bench, fc, mc_analyse as ma
S = bench.load_set('w1024_d16'); W = bench.weights(S, 0)
tr = []; fc.run(W, trace=tr, slices=2)
D = np.load('data/mc1024_mlp0_N65536_a.npz')
def cen(h, l):
    d = D[f'{h}_s1t'][l]; E2 = D[f'{h}_s2t'][l].astype(float); E21 = D[f'{h}_s21a'][l].astype(float); E3 = D[f'{h}_s3a'][l]
    e2 = np.diag(E2)
    K = E21 - 2 * E2 * d[:, None] - e2[:, None] * d[None, :] + 2 * (d * d)[:, None] * d[None, :]
    k3 = E3 - 3 * d * e2 + 2 * d ** 3
    np.fill_diagonal(K, 0)
    return K, k3
def eps(M, A, B):
    return np.sqrt(max(np.sum((M - A) * (M - B)), 0) / max(np.sum(A * B), 1e-30))
for l in range(1, 8):
    A, a3 = cen('A', l); B, b3 = cen('B', l)
    t = tr[l]
    print(l, 'K21: model', f"{eps(t['sl21'], A, B):.3f}", 'gauss', f"{eps(t['sl21G'], A, B):.3f}",
          ' k3: model', f"{eps(t['sl3'], a3, b3):.3f}", 'gauss', f"{eps(t['sl3G'], a3, b3):.3f}",
          ' |K21|', f"{np.sqrt(max(np.sum(A*B),0)):.3f}", flush=True)
