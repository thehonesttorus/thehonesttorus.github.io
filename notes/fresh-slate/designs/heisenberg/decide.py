"""Deciding experiment: computed first-order HD + TRUE joint kappa_4 slices (MC), full second-order injection."""
import numpy as np, sys, json
sys.path.insert(0, '/root/hdw'); sys.path.insert(0, __file__.rsplit('/', 1)[0])
from oracle2 import passes
from hd import hd
n = int(sys.argv[1]); N = int(float(sys.argv[2])); L = 16
for seed in range(int(sys.argv[3]), int(sys.argv[4])):
    rng = np.random.default_rng(500 + seed); Ws = rng.normal(size=(L, n, n)) * np.sqrt(2 / n)
    p1 = passes(Ws, N, seed + 50); st = passes(Ws, N, seed + 60, mean=p1['m']); tr = p1['a']
    k4 = [None]
    for l in range(1, L):
        c2 = st['c2'][l]; v = np.diag(c2)
        k4.append((st['d4'][l] - 3 * v * v, st['s22'][l] - np.outer(v, v) - 2 * c2 * c2, st['s31'][l] - 3 * v[:, None] * c2))
    h1 = hd(Ws, A=None, Kt=3, diagrams="star")
    h2 = hd(Ws, A=None, Kt=3, diagrams="star", k4ext=k4)
    e1 = ((h1 - tr) ** 2).mean(1); e2 = ((h2 - tr) ** 2).mean(1)
    print(json.dumps(dict(n=n, seed=seed, hd1=e1[-1], hd1_trueK4=e2[-1], gain=e1[-1] / e2[-1], noise=0.15 / N)), flush=True)
