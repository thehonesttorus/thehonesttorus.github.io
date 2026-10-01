"""Oracle injection at depth 16: replace chosen parts of the Bethe state by Monte Carlo truth at every
layer (pre-ReLU z: m, C, K = kappa(z_a,z_a,z_b) incl. kappa3, kappa4) and measure the final MSE."""
import sys, numpy as np
sys.path.insert(0, '../../bench'); sys.path.insert(0, '.')
import bench, bethe


def mc_state(W, N, batch=100_000, seed=7):
    rng = np.random.default_rng(seed)
    L, n, _ = W.shape
    s1 = np.zeros((L, n)); s2 = np.zeros((L, n, n)); s21 = np.zeros((L, n, n)); s4 = np.zeros((L, n)); s3 = np.zeros((L, n))
    done = 0
    while done < N:
        a = rng.standard_normal((batch, n))
        for l in range(L):
            z = a @ W[l]; a = np.maximum(z, 0)
            s1[l] += z.sum(0); s2[l] += z.T @ z; s21[l] += (z * z).T @ z; s3[l] += (z ** 3).sum(0); s4[l] += (z ** 4).sum(0)
        done += batch
    out = []
    for l in range(L):
        m = s1[l] / N; E2 = s2[l] / N; E21 = s21[l] / N
        C = E2 - np.outer(m, m)
        # kappa(z_a,z_a,z_b) = E[a^2 b] - E[a^2]E[b] - 2E[a](E[ab]-E[a]E[b])
        K = E21 - np.outer(np.diag(E2), m) - 2 * m[:, None] * C
        v = np.diag(C); Ez2 = np.diag(E2); Ez3 = s3[l] / N
        k4 = s4[l] / N - 4 * Ez3 * m - 3 * Ez2 ** 2 + 12 * Ez2 * m * m - 6 * m ** 4
        out.append(dict(m=m, C=C, K=K, k4=k4))
    return out


def run(Ws, truth, inject, old=1):
    Ls, n, _ = Ws.shape
    W = Ws[0]
    m = np.zeros(n); C = W.T @ W; K = np.zeros((n, n)); k4 = np.zeros(n)
    out = []; prev = None
    for l in range(Ls):
        t = truth[l]
        if 'm' in inject: m = t['m']
        if 'C' in inject: C = t['C']
        if 'v' in inject: C = C.copy(); np.fill_diagonal(C, np.diag(t['C']))
        if 'K' in inject: K = t['K']
        if 'k3' in inject: K = K.copy(); np.fill_diagonal(K, np.diag(t['K']))
        if 'k4' in inject: k4 = t['k4']
        mu, Ca, Ka, k3a, k4a, c, L0, Lm1 = bethe.relu_map_edges(m, C, K, k4)
        out.append(mu)
        if l + 1 == Ls:
            break
        Wn = Ws[l + 1]
        M = c * c * (L0 ** 2)[:, None] * Lm1[None, :]
        spec = (k3a, Ka - M, c, L0, Lm1)
        Kn = bethe.contract(*spec, Wn, Wn)
        if old and prev is not None:
            pspec, Wl = prev
            P = Wl @ (L0[:, None] * Wn)
            Kz = K.copy(); k3z = np.diag(K).copy(); np.fill_diagonal(Kz, 0.0)
            Kn = Kn + bethe.contract(*pspec, P, P) - bethe.contract(k3z * L0 ** 3, Kz * (L0 ** 2)[:, None] * L0[None, :],
                                                                  None, L0, Lm1, Wn, Wn)
        prev = (spec, Wn)
        m = mu @ Wn; C = Wn.T @ Ca @ Wn; K = Kn; k4 = (Wn ** 4).T @ k4a
    return np.array(out)


if __name__ == '__main__':
    name, N, mlps = sys.argv[1], float(sys.argv[2]), [int(x) for x in sys.argv[3].split(',')]
    S = bench.load_set(name)
    sets = [(), ('k4',), ('k3',), ('k3', 'k4'), ('v',), ('v', 'k3', 'k4'), ('K',), ('K', 'k4'), ('C',), ('C', 'k4'),
            ('C', 'K'), ('C', 'K', 'k4'), ('m', 'C', 'K', 'k4')]
    print(name, 'N', N, 'inject:', ' | '.join('+'.join(s) or 'none' for s in sets), flush=True)
    for i in mlps:
        W = bench.weights(S, i).astype(np.float64)
        tr = mc_state(W, int(N))
        T = S['means'][i]
        row = []
        for s in sets:
            e = run(W, tr, s)
            row.append(((e[-1] - T[-1]) ** 2).mean())
        print(f"mlp {i}: " + ' '.join(f'{x:.2e}' for x in row), flush=True)
