"""Rank-versus-age law of the gated propagators Z_s(t) = W_{s+1} D_{Phi_{s+1}} ... W_t (Gaussian-closure gates), the
object whose frozen-frame resolution sets the carrier cost (team D: k(a) = 2n/a lossless at 1024 x 16).
For every (s, t): participation ratio PR, and k_q = number of singular directions holding a fraction q of ||Z||_F^2.
Fit k(a) ~ a^{-p} over ages a = t - s; p = 1 is Dixmier-critical (cost ~ L log L), p < 1 gives L^{2-p}, p > 1 gives L.
usage: python g3_ranks.py SET MLP"""
import json, sys, os, numpy as np
from scipy.special import ndtr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '../../fresh-slate/bench')); sys.path.insert(0, os.path.join(HERE, '../../fresh-slate/breakthrough/region'))
import bench, gclose as g
name, i = sys.argv[1], int(sys.argv[2])
S = bench.load_set(name); W = bench.weights(S, i).astype(np.float64); L, n, _ = W.shape
_, states = g.run(W, 'exact')
Phi = []
for l in range(L):
    if l == 0: mu = np.zeros(n); Sg = W[0].T @ W[0]
    else: m, C = states[l - 1]; mu = m @ W[l]; Sg = W[l].T @ C @ W[l]
    Phi.append(ndtr(mu / np.sqrt(np.diag(Sg))))
rows = []
for s in range(L - 1):
    Z = W[s + 1].copy()                       # Z_s(s+1)
    for t in range(s + 1, L):
        if t > s + 1:
            Z = (Z * Phi[t - 1][None, :]) @ W[t]
        sv2 = np.linalg.svd(Z, compute_uv=False) ** 2; tot = sv2.sum(); c = np.cumsum(sv2) / tot
        rows.append(dict(s=s, t=t, a=t - s, PR=float(tot ** 2 / (sv2 ** 2).sum()),
                         k90=int(np.searchsorted(c, 0.90) + 1), k99=int(np.searchsorted(c, 0.99) + 1), top=float(sv2[0] / tot)))
json.dump(rows, open(os.path.join(HERE, f'results/g3ranks_{name}_{i}.json'), 'w'))
A = np.array([r['a'] for r in rows])
print(f'{name} mlp {i}: n={n} L={L}')
print(' a   PR/n   k90/n  k99/n  top-share   (means over s)   n/(2(a+1))/n')
for a in sorted(set(A)):
    sel = [r for r in rows if r['a'] == a]
    if a <= 6 or a % 4 == 0 or a == max(A):
        print(f"{a:2d}  {np.mean([r['PR'] for r in sel])/n:.3f}  {np.mean([r['k90'] for r in sel])/n:.3f}  {np.mean([r['k99'] for r in sel])/n:.3f}  {np.mean([r['top'] for r in sel]):.3f}      {1/(2*(a+1)):.3f}")
for key in ('PR', 'k90', 'k99'):
    for lo, hi in ((2, 8), (8, L - 1)):
        sel = [(r['a'], r[key]) for r in rows if lo <= r['a'] <= hi]
        if len(sel) > 3:
            x = np.log([a for a, _ in sel]); y = np.log([v for _, v in sel]); p = -np.polyfit(x, y, 1)[0]
            print(f'  fit {key} ~ a^-p over ages {lo}-{hi}: p = {p:.2f}')
