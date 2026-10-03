# Depth-6 cells with more networks: per-network boundary mass, per-layer |.|-mass, mean gradient norm^2.
import numpy as np
from multiprocessing import Pool
from edge import net, chunk_stats
d, CH, J = 8, 50_000, 8
def work(args):
    n, L, seed = args
    Ws = net(d, n, L, np.random.default_rng(seed))
    tot = None
    for c in range(4):
        X = np.random.default_rng(10_000*seed + c).standard_normal((CH, d))
        o = chunk_stats(Ws, X, list(range(J)), (0.04, 0.02))
        A = np.array([2*o[j][3][1] - o[j][3][0] for j in range(J)])/CH
        m = np.array([o[j][0] for j in range(J)])/CH
        tot = (A, m) if tot is None else (tot[0] + A, tot[1] + m)
    A, m = tot[0]/4, tot[1]/4
    # gradient-norm check: E |grad z_{L-1,i}|^2 averaged over neurons
    X = np.random.default_rng(seed).standard_normal((20_000, d))
    H = X; Jp = np.broadcast_to(np.eye(d), (len(X), d, d))
    for W in Ws[:-1]:
        Z = H @ W.T; Jz = np.einsum('ab,nbc->nac', W, Jp); H = np.maximum(Z, 0); Jp = (Z > 0)[:, :, None]*Jz
    g2 = (Jz**2).sum(2).mean()
    return A.sum(1).mean(), A[:, :-1].mean(), m.mean(), g2
if __name__ == '__main__':
    with Pool(4) as p:
        for n in [32, 64]:
            L = 6
            R = np.array(p.map(work, [(n, L, s) for s in range(100, 108)]))
            print(f"n={n} L={L}: per-network ||LapF|| {np.round(R[:,0],2)} mean {R[:,0].mean():.2f} (heur {1/np.sqrt(np.pi)+(L-1)*np.sqrt(n)/np.pi:.2f});"
                  f" per-layer |.| mean {R[:,1].mean():.3f} (heur {np.sqrt(n)/np.pi:.3f}); E F mean {R[:,2].mean():.3f}; |grad z_(L-1)|^2 per net {np.round(R[:,3],2)}", flush=True)
