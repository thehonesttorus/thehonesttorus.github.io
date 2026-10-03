# Larger run of the edge formula with chunk-level standard errors (4 processes).
import numpy as np, sys
from multiprocessing import Pool
from edge import net, chunk_stats
d, n, L = 8, 32, 4
js = [0, 1, 2, 3, 4, 5]
epss = (0.04, 0.02, 0.01)
CH = 100_000
def work(c):
    rng = np.random.default_rng(1); Ws = net(d, n, L, rng)       # same network as edge.py
    X = np.random.default_rng(1000 + c).standard_normal((CH, d))
    o = chunk_stats(Ws, X, js, epss)
    return {j: (o[j][0]/CH, o[j][2].sum(1)/CH, o[j][3].sum(1)/CH, o[j][2]/CH) for j in js}
if __name__ == '__main__':
    nch = int(sys.argv[1]) if len(sys.argv) > 1 else 120
    with Pool(4) as p: R = p.map(work, range(nch))
    for j in js:
        mc = np.array([r[j][0] for r in R]); S = np.array([r[j][1] for r in R]); A = np.array([r[j][2] for r in R])
        Pl = np.array([r[j][3] for r in R])                          # (chunks, eps, L)
        rich_a = 2*S[:, 1] - S[:, 0]; rich_b = 2*S[:, 2] - S[:, 1]
        f = lambda v: f"{v.mean():.4f}+-{v.std()/np.sqrt(len(v)):.4f}"
        own = 2*Pl[:, 1, L-1] - Pl[:, 0, L-1]
        deeper = (2*Pl[:, 1, :L-1] - Pl[:, 0, :L-1]).sum(1)
        tv = (2*A[:, 1] - A[:, 0]).mean()
        print(f"j={j}: MC {f(mc)} | edge (.04/.02) {f(rich_a)} | edge (.02/.01) {f(rich_b)} | own kink {f(own)} earlier layers {f(deeper)} | ||Lap F|| {tv:.3f} rho {mc.mean()/tv:.3f}")
