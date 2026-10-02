"""Linear-response fit of birth-coefficient changes against MC D21 truth, per layer.
D(c) ~ D0 + sum_k dc_k Delta_k with Delta_k from global perturbation runs. Off-diagonal entries; split-half truth.
Prints per layer: eps(base), optimal dc_k (in the knob's units), predicted eps after the fit."""
import numpy as np
mc = np.load('mc_d21_m0_n19.npz'); A = mc['D21_A'].astype(np.float64); Bm = mc['D21_B'].astype(np.float64)
n = A.shape[1]; off = ~np.eye(n, dtype=bool)
base = np.load('rec2_F1final_m0.npz')                      # FB 1.5, LAM 0.80, FEED 1.0
runs = {'FB': (np.load('rec2_F1final_fb1_m0.npz'), 1.0 - 1.5),
        'LAM': (np.load('rec2_F1_lam95_m0.npz'), 0.95 - 0.80),
        'FEED': (np.load('rec2_F1_feed13_m0.npz'), 1.3 - 1.0)}
print(' l   eps0    ' + '  '.join('%-8s' % k for k in runs) + '  eps_fit   (dc in knob units: FB scale, LAM scale, FEED scale)')
for l in range(2, 15):
    D0 = base[f'D21_{l}'].astype(np.float64)[off]; a = A[l][off]; b = Bm[l][off]; t = (a + b) / 2
    sig = (a * b).sum()
    X = np.stack([(r[f'D21_{l}'].astype(np.float64)[off] - D0) / h for r, h in runs.values()], 1)
    G = X.T @ X; g = X.T @ (t - D0)
    dc = np.linalg.solve(G + 1e-12 * np.trace(G) * np.eye(len(runs)), g)
    Df = D0 + X @ dc
    e0 = np.sqrt(max(((D0 - a) * (D0 - b)).sum(), 0) / sig); ef = np.sqrt(max(((Df - a) * (Df - b)).sum(), 0) / sig)
    print('%2d  %6.4f  ' % (l, e0) + '  '.join('%+8.3f' % x for x in dc) + '  %6.4f' % ef)
