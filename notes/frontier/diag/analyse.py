"""Per-layer error split (bias^2 vs rest) and leave-one-MLP-out affine calibration of the final layer.
usage: python3 analyse.py pred.npz"""
import sys
import numpy as np
d = np.load(sys.argv[1]); P, T, NZ = d['pred'], d['truth'], d['noise']
m, L, n = P.shape
E = P - T
print('final-layer raw per MLP (x1e-8):', ' '.join('%.3f' % ((((E[i, -1]) ** 2).mean() - NZ[i]) * 1e8) for i in range(m)))
print('layer  bias^2/MSE  (mean over MLPs; MSE incl. truth noise)')
for l in range(L):
    b2 = np.mean([E[i, l].mean() ** 2 for i in range(m)]); mse = np.mean([(E[i, l] ** 2).mean() for i in range(m)])
    print('%2d  bias^2 %.2e  mse %.2e  ratio %.3f  mean bias %+.2e' % (l, b2, mse, b2 / mse, np.mean([E[i, l].mean() for i in range(m)])))
# affine calibration of the last layer, leave one MLP out
raw0, raw1, raw2 = [], [], []
for i in range(m):
    tr = [j for j in range(m) if j != i]
    X = np.concatenate([np.stack([np.ones(n), P[j, -1]], 1) for j in tr]); y = np.concatenate([-E[j, -1] for j in tr])
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    corr = beta[0] + beta[1] * P[i, -1]
    raw0.append((E[i, -1] ** 2).mean() - NZ[i]); raw1.append(((E[i, -1] + corr) ** 2).mean() - NZ[i])
    c2 = -E[i, -1].mean() * 0 + np.mean([E[j, -1].mean() for j in tr])   # pooled constant-bias removal
    raw2.append(((E[i, -1] - c2) ** 2).mean() - NZ[i])
print('LOO affine calib: raw %.4e -> %.4e (%+.1f%%); constant-bias removal -> %.4e (%+.1f%%)' % (
    np.mean(raw0), np.mean(raw1), 100 * (np.mean(raw1) / np.mean(raw0) - 1), np.mean(raw2), 100 * (np.mean(raw2) / np.mean(raw0) - 1)))
