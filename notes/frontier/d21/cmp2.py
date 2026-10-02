"""Per-layer errors vs split-half MC truth: D21 off-diagonal, D3 = diag kappa3, var; unbiased squared-norm estimators.
usage: python3 cmp2.py MC.npz REC.npz [REC2 ...]"""
import sys
import numpy as np
mc = np.load(sys.argv[1]); A = mc['D21_A'].astype(np.float64); Bm = mc['D21_B'].astype(np.float64)
vA, vB = mc['var_A'], mc['var_B']; L, n, _ = A.shape
off = ~np.eye(n, dtype=bool); di = np.arange(n)
def rel(est, a, b):
    sig = (a * b).sum(); err = ((est - a) * (est - b)).sum()
    s = (est * (a + b) / 2).sum() / (est * est).sum()
    errs = ((s * est - a) * (s * est - b)).sum()
    return np.sqrt(max(err, 0) / sig), s, np.sqrt(max(errs, 0) / sig)
for f in sys.argv[2:]:
    r = np.load(f); print('==', f)
    print(' l | D21 off: eps  scale eps|s | D3: eps  scale | var: eps(x1e3)  scale-1(x1e3)')
    for l in range(1, L):
        if f'D21_{l}' not in r: continue
        D = r[f'D21_{l}'].astype(np.float64)
        e1, s1, e1s = rel(D[off], A[l][off], Bm[l][off])
        d3 = r[f'D3_{l}'].astype(np.float64) if f'D3_{l}' in r else np.diag(D)
        e2, s2, _ = rel(d3, A[l][di, di], Bm[l][di, di])
        v = r[f'var_{l}'].astype(np.float64) if f'var_{l}' in r else None
        if v is not None:
            e3, s3, _ = rel(v, vA[l], vB[l])
        else:
            e3, s3 = np.nan, np.nan
        print('%2d | %6.3f %6.3f %6.3f | %6.3f %6.3f | %7.3f %7.3f' % (l, e1, s1, e1s, e2, s2, 1e3 * e3, 1e3 * (s3 - 1)))
