"""Proposition E1 (wedge calculus) against linear-response probes at other shapes (results/lr_<set>_mlp<i>.json,
from lr_probe_any.py) and region's n = 1024 x 16 probes.  No fitted parameter."""
import json, glob, os, numpy as np

def wedge(n, L):
    f = lambda r: (np.sqrt(1 - r * r) + (np.pi - np.arccos(r)) * r) / np.pi
    rho = [0.0]
    for _ in range(L - 1):
        rho.append(f(rho[-1]))
    rho = np.array(rho); lam = 0.5 + np.arcsin(rho) / np.pi
    Km = np.array([np.prod(lam[s + 1:L]) for s in range(L)])
    q = (1 / (2 * np.pi)) / np.sqrt(1 + 2 * rho / (1 - rho)) / (4 * 2 * (1 - rho))
    Ko = np.zeros(L)
    for l in range(L - 1):
        Ko[l] = (8 / n ** 2) * sum(np.prod(lam[l + 1:k] ** 2) * q[k] * Km[k] for k in range(l + 1, L))
    return Km, Ko

def compare(tag, files):
    M = [json.load(open(f)) for f in files]
    L = len(M[0]['mean']); n = M[0].get('n', 1024)
    Km, Ko = wedge(n, L)
    lm = np.array([[np.log(m['mean'][l] / Km[l]) for m in M] for l in range(L - 1)])
    lo = np.array([[np.log(m['off'][l] / Ko[l]) for m in M] for l in range(L - 1)])
    print(f'{tag}: n={n} L={L} MLPs={len(M)}  log(meas/pred) K_mean: mean {lm.mean():+.2f} rms {np.sqrt((lm**2).mean()):.2f}'
          f' | K_off: mean {lo.mean():+.2f} rms {np.sqrt((lo**2).mean()):.2f}')
    for l in list(range(0, L - 1, max(1, (L - 1) // 8))) + [L - 2]:
        print(f'   l={l:2d}  Kmean pred {Km[l]:.3f} meas ' + ' '.join(f"{m['mean'][l]:.3f}" for m in M)
              + f'   Koff pred {Ko[l]:.2e} meas ' + ' '.join(f"{m['off'][l]:.2e}" for m in M))

reg = '../../fresh-slate/breakthrough/region/results/'
compare('w1024_d16 (region)', [reg + f'lr_mlp{i}.json' for i in range(3)])
for name in ['w512_d16', 'w256_d32', 'w128_d16']:
    fs = sorted(glob.glob(f'results/lr_{name}_mlp*.json'))
    if fs:
        compare(name, fs)
