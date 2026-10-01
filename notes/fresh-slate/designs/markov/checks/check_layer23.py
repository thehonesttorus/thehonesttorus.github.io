"""High-precision MC of z_2, z_3 cumulants (64 wide, first 3 layers only) vs MKV."""
import sys, numpy as np
sys.path.insert(0, '.')
from mkv import mkv
from bake import weights
n = int(sys.argv[1]); N = int(float(sys.argv[2]))
W = weights(n, 16, 1)[:3].astype(np.float64)
rng = np.random.default_rng(5)
S = np.zeros((2, 4, n)); cnt = 0
while cnt < N:
    x = rng.standard_normal((200000, n))
    z2 = np.maximum(x @ W[0], 0) @ W[1]; z3 = np.maximum(z2, 0) @ W[2]
    for i, z in enumerate([z2, z3]):
        p = np.ones_like(z)
        for k in range(4):
            p = p * z; S[i, k] += p.sum(0)
    cnt += len(x)
m = S / cnt
_, dg = mkv(W, w=4, var21='full', return_all=True)
for i in range(2):
    m1, m2, m3, m4 = m[i]
    v = m2 - m1**2; k3 = m3 - 3*m2*m1 + 2*m1**3; k4 = m4 - 4*m3*m1 - 3*m2**2 + 12*m2*m1**2 - 6*m1**4
    d = dg[i + 1]
    r = lambda a, b: np.sqrt(np.mean((a - b)**2)) / np.sqrt(np.mean(b**2))
    print(f"layer {i+2}: var {r(d['v'], v):.2e}  k3 {r(d['k3'], k3):.2e}  k4 {r(d['k4'], k4):.2e}  corr(k3) {np.corrcoef(d['k3'], k3)[0,1]:.4f}  slope {np.dot(d['k3'],k3)/np.dot(k3,k3):.3f}")
