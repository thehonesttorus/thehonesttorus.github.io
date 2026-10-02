"""Per-layer D21 error of an estimator against split-half MC truth (noise-free estimators of the squared norms).
usage: python3 cmp_d21.py MC.npz REC.npz [REC2.npz ...]"""
import sys
import numpy as np
mc = np.load(sys.argv[1]); A = mc['D21_A'].astype(np.float64); Bm = mc['D21_B'].astype(np.float64); L, n, _ = A.shape
off = ~np.eye(n, dtype=bool)
def stats(De, l, orient):
    a, b = A[l], Bm[l]
    D = De if orient == 0 else De.T
    sig = (a[off] * b[off]).sum()                       # ||D_true||^2 (off-diagonal), unbiased
    err = ((D - a)[off] * (D - b)[off]).sum()           # ||D - D_true||^2, unbiased
    s = (D[off] * (a + b)[off] / 2).sum() / (D[off] ** 2).sum()   # best scale
    errs = ((s * D - a)[off] * (s * D - b)[off]).sum()
    return sig, err, s, errs
for f in sys.argv[2:]:
    r = np.load(f)
    print('==', f)
    print(' l   eps(off)   best-scale  eps-after-scale   ||D_true||_F   orient')
    for l in range(1, L):
        k = f'D21_{l}'
        if k not in r: continue
        De = r[k].astype(np.float64)
        best = min((stats(De, l, o) + (o,) for o in (0, 1)), key=lambda t: t[1])
        sig, err, s, errs, o = best
        print('%2d   %7.3f    %7.3f      %7.3f          %.3e     %d' % (l, np.sqrt(max(err, 0) / sig), s, np.sqrt(max(errs, 0) / sig), np.sqrt(max(sig, 0)), o))
