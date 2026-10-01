"""Split FC's D21(z_{l+1}) into (old) the gated transport of kappa3(z_l) and (birth) the layer-l source, and check
each against Monte Carlo at n = 1024: MC old = E[x_a^2 x_b] (central), x = (P_l o (z_l - mu_l)) W_{l+1} with the model's
gates P_l; MC birth = MC D21(z_{l+1}) - MC old.  Two halves, cross-product eps.
usage: python check_b0.py MLP N_HALF LMAX"""
import sys
import numpy as np
sys.path.insert(0, '../../bench')
import bench, fc  # noqa: E402
from scipy.special import ndtr


def cen21(E21, E2, d):
    e2 = np.diag(E2)
    return E21 - 2 * E2 * d[:, None] - e2[:, None] * d[None, :] + 2 * (d * d)[:, None] * d[None, :]


def main(i, nh, lmax):
    S = bench.load_set('w1024_d16'); W = bench.weights(S, i); T = S['means'][i]
    L, n, _ = W.shape
    # FC trace with per-source split at layers 1..lmax
    tr = []
    split = {}
    orig = fc.run
    fc_kw = dict(slices=2)
    # rerun FC manually to capture the split: windowed runs give D21 from sources born at l only (window=1)
    p_all = fc.run(W, trace=tr, **fc_kw)
    tr1 = []
    fc.run(W, trace=tr1, window=1, **fc_kw)
    Ps = []
    for l in range(L):
        mu, Sm = tr[l]['mu'], tr[l]['S']
        Ps.append(ndtr(mu / np.sqrt(np.diag(Sm))))
    mu_ref = np.zeros((L, n))
    for l in range(1, L):
        mu_ref[l] = T[l - 1] @ W[l].astype(np.float64)
    acc = {}
    for h in (0, 1):
        rng = np.random.default_rng([555 + i, h])
        done = 0
        while done < nh:
            m = 4096
            x0 = rng.standard_normal((m, n), dtype=np.float32)
            a = x0
            for l in range(lmax + 1):
                z = a @ W[l]
                y = z - mu_ref[l].astype(np.float32)
                a = np.maximum(z, 0)
                if l < lmax:
                    xl = (y * Ps[l].astype(np.float32)) @ W[l + 1]   # x for target l+1
                    zn = a @ W[l + 1]; yn = zn - mu_ref[l + 1].astype(np.float32)
                    for key, v in (('x', xl), ('y', yn)):
                        k = (h, l, key)
                        st = acc.setdefault(k, [0, 0, 0])
                        st[0] = st[0] + v.sum(0, dtype=np.float64)
                        st[1] = st[1] + (v.T @ v).astype(np.float64)
                        st[2] = st[2] + ((v * v).T @ v).astype(np.float64)
            done += m
    print('layer(target)  eps_total  eps_old  eps_birth  share_old(model)  |old|/|D21|(MC)')
    for l in range(lmax):
        res = {}
        for key in ('x', 'y'):
            for h in (0, 1):
                s1, s2, s21 = [v / nh for v in acc[(h, l, key)]]
                res[(key, h)] = cen21(s21, s2, s1)
        Dt = tr[l + 1]['D21']; Db = tr1[l + 1]['D21']; Do = Dt - Db
        def eps(M, A, B):
            return np.sqrt(max(np.sum((M - A) * (M - B)), 0) / np.sum(A * B))
        oA, oB = res[('x', 0)], res[('x', 1)]
        tA, tB = res[('y', 0)], res[('y', 1)]
        print(l + 1, f"{eps(Dt, tA, tB):.3f} {eps(Do, oA, oB):.3f} {eps(Db, tA - oA, tB - oB):.3f} "
              f"{np.linalg.norm(Do) / np.linalg.norm(Dt):.3f} {np.sqrt(np.sum(oA * oB) / np.sum(tA * tB)):.3f}", flush=True)


if __name__ == '__main__':
    main(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]))
