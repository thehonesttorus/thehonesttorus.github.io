"""Teacher forcing: apply the Bethe ReLU map to the MC-true z-state of each layer and compare its
output with the MC-true a-state.  Isolates the closure error of each piece."""
import sys, numpy as np
sys.path.insert(0, '.')
from bethe import *


def cum3(E1, E2, E3):
    return (E3 - np.einsum('ij,k->ijk', E2, E1) - np.einsum('ik,j->ijk', E2, E1)
            - np.einsum('jk,i->ijk', E2, E1) + 2 * np.einsum('i,j,k->ijk', E1, E1, E1))


def mc_stats(Ws, N, batch=50_000, seed=1):
    rng = np.random.default_rng(seed)
    Ls, n, _ = Ws.shape
    acc = [dict() for _ in range(Ls)]
    done = 0
    while done < N:
        a = rng.standard_normal((batch, n))
        for l in range(Ls):
            z = a @ Ws[l]; a = np.maximum(z, 0)
            d = acc[l]
            for nm, x in (('z', z), ('a', a)):
                X2 = (x[:, :, None] * x[:, None, :]).reshape(len(x), -1)
                for k, val in (('1', x.sum(0)), ('2', x.T @ x),
                               ('3', (X2.T @ x).reshape(n, n, n))):
                    d[nm + k] = d.get(nm + k, 0) + val
            d['z4'] = d.get('z4', 0) + (z ** 4).sum(0)
        done += batch
    out = []
    for d in acc:
        st = {}
        for nm in 'za':
            E1 = d[nm + '1'] / N; E2 = d[nm + '2'] / N; E3 = d[nm + '3'] / N
            st[nm] = (E1, E2 - np.outer(E1, E1), cum3(E1, E2, E3))
        m = st['z'][0]; v = np.diag(st['z'][1]); k3 = np.einsum('aaa->a', st['z'][2])
        Ez2 = v + m * m; Ez3 = k3 + 3 * m * v + m ** 3
        st['k4'] = d['z4'] / N - 4 * Ez3 * m - 3 * Ez2 ** 2 + 12 * Ez2 * m * m - 6 * m ** 4
        out.append(st)
    return out


if __name__ == '__main__':
    n, L, N = int(sys.argv[1]), int(sys.argv[2]), int(float(sys.argv[3]))
    rng = np.random.default_rng(int(sys.argv[4]) if len(sys.argv) > 4 else 0)
    Ws = rng.standard_normal((L, n, n)) * np.sqrt(2 / n)
    st = mc_stats(Ws, N)
    idx = np.arange(n); off = ~np.eye(n, dtype=bool)
    dist = np.ones((n, n, n), bool); dist[idx, idx, :] = 0; dist[idx, :, idx] = 0; dist[:, idx, idx] = 0
    rms = lambda x: np.sqrt(np.mean(x ** 2))
    print(f'n={n} N={N:.0e}; layer | mu err (gauss, +k3, +k3k4) | C^a off rel | K^a off rel | '
          'T^a rel (hubs+old, hubs only) | rms c, K_z, T_z')
    for l in range(L):
        mz, Cz, Sz = st[l]['z']; ma, Ca_t, Sa_t = st[l]['a']; k4 = st[l]['k4']
        r = []
        for k3s, k4s in ((0, 0), (1, 0), (1, 1)):
            mu, _, _, _ = relu_map_full(mz, Cz, Sz * k3s, k4 * k4s)
            r.append(rms(mu - ma))
        mu, Ca, Sa, _ = relu_map_full(mz, Cz, Sz, k4)
        _, _, Sa0, _ = relu_map_full(mz, Cz, Sz, k4, keep_old=False)
        Kt = Sa_t[idx, idx, :]; K = Sa[idx, idx, :]
        print(f"{l:2d} | {r[0]:.1e} {r[1]:.1e} {r[2]:.1e} | {rms((Ca-Ca_t)[off])/rms(Ca_t[off]):.3f} | "
              f"{rms((K-Kt)[off])/rms(Kt[off]):.3f} | {rms((Sa-Sa_t)[dist])/rms(Sa_t[dist]):.3f} "
              f"{rms((Sa0-Sa_t)[dist])/rms(Sa_t[dist]):.3f} | {rms(Cz[off]):.3f} "
              f"{rms(Sz[idx,idx,:][off]):.4f} {rms(Sz[dist]):.5f}", flush=True)
