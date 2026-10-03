# Boundary mass ||Lap F||_gamma = sum_e gamma_{d-1}(e)|kappa_e| versus depth L and width n, against the
# mean-field heuristic 1/sqrt(pi) + (L-1) sqrt(n)/pi; cancellation index rho = E F / ||Lap F||.
# Also: own-kink term vs the Gaussian-width proxy sqrt(v) phi(m/sqrt v), m = E z_Lj, v = Var z_Lj.
import numpy as np, sys
from math import erf
from multiprocessing import Pool
from edge import net, chunk_stats
d, CH, J = 8, 50_000, 8
def work(args):
    n, L, seed, c = args
    Ws = net(d, n, L, np.random.default_rng(seed))
    X = np.random.default_rng(10_000*seed + c).standard_normal((CH, d))
    o = chunk_stats(Ws, X, list(range(J)), (0.04, 0.02))
    mean = np.array([o[j][0] for j in range(J)])/CH
    A = np.array([2*o[j][3][1] - o[j][3][0] for j in range(J)])/CH     # (J, L)
    S = np.array([2*o[j][2][1] - o[j][2][0] for j in range(J)])/CH
    z1 = np.array([o[j][5] for j in range(J)])/CH; z2 = np.array([o[j][6] for j in range(J)])/CH
    return mean, A, S, z1, z2
phi = lambda t: np.exp(-t*t/2)/np.sqrt(2*np.pi)
Phi = np.vectorize(lambda t: 0.5*(1+erf(t/np.sqrt(2))))
if __name__ == '__main__':
    nch = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    with Pool(4) as p:
        for n in [16, 32, 64]:
            for L in [2, 3, 4, 6]:
                R = p.map(work, [(n, L, s, c) for s in (11, 12) for c in range(nch)])
                res = []
                for s in range(2):
                    Rs = R[s*nch:(s+1)*nch]
                    mean, A, S, z1, z2 = [np.mean([r[k] for r in Rs], axis=0) for k in range(5)]
                    res.append((mean, A, S, z1, z2 - z1**2))
                mean = np.concatenate([r[0] for r in res]); A = np.concatenate([r[1] for r in res]); S = np.concatenate([r[2] for r in res])
                m = np.concatenate([r[3] for r in res]); v = np.concatenate([r[4] for r in res])
                tv = A.sum(1); own = S[:, -1]; rest = S[:, :-1].sum(1)
                width_proxy = np.sqrt(v)*phi(m/np.sqrt(v)); shift_proxy = m*Phi(m/np.sqrt(v))
                pred = 1/np.sqrt(np.pi) + (L-1)*np.sqrt(n)/np.pi
                print(f"n={n:3d} L={L}: E F {mean.mean():.3f} | ||LapF|| {tv.mean():.3f} (heur {pred:.3f}) | rho med {np.median(mean/tv):.4f} |"
                      f" own {own.mean():.3f} vs width proxy {width_proxy.mean():.3f} (corr {np.corrcoef(own, width_proxy)[0,1]:.2f}) |"
                      f" earlier {rest.mean():.3f} vs shift proxy {shift_proxy.mean():.3f} (corr {np.corrcoef(rest, shift_proxy)[0,1]:.2f}) |"
                      f" per-layer |.| {np.round(A.mean(0),3)}", flush=True)
