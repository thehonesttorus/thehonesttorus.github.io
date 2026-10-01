"""kappa_3 slice accuracy of the pull-back HD with mean-gate transport vs with the exact coincident (re-birth) response."""
import numpy as np, sys
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from oracle2 import passes
from hdpull import hd_pull
n = int(sys.argv[1]) if len(sys.argv) > 1 else 64; L = 16
rng = np.random.default_rng(500); Ws = rng.normal(size=(L, n, n)) * np.sqrt(2 / n)
p1 = passes(Ws, int(float(sys.argv[2])) if len(sys.argv) > 2 else 3_000_000, 50); st = passes(Ws, int(float(sys.argv[2])) if len(sys.argv) > 2 else 3_000_000, 60, mean=p1['m'])
off = ~np.eye(n, dtype=bool)
for rb in [False, True]:
    tr = []; hd_pull(Ws, rebirth=rb, trace=tr)
    for l, D, S in tr:
        if l in (1, 2, 4, 8, 12, 15):
            Dt = st['d3'][l]; St = st['s21'][l]
            print("rebirth" if rb else "meangate", l, "D rel %.3f | S rel %.3f" % (np.linalg.norm(D - Dt) / np.linalg.norm(Dt), np.linalg.norm((S - St)[off]) / np.linalg.norm(St[off])), flush=True)
