"""Isolate the all-distinct birth at layer lb (target lb+1) and compare FC's leading-Wick tensor with MC.
MC birth_ad = D21(z_{lb+1}) - T[Phi^3 kappa3(z_lb)] - T[slices(kappa3(a_lb)) - Phi^3 slices(kappa3(z_lb))]
where T is the transport by W_{lb+1} of a slice-supported tensor (exact formula), slices from the a-atlas.
usage: python check_birth.py LB N_HALF"""
import sys
import numpy as np
from scipy.special import ndtr
sys.path.insert(0, '../../bench')
import bench, fc  # noqa: E402


def cen21(E21, E2, d):
    e2 = np.diag(E2)
    return E21 - 2 * E2 * d[:, None] - e2[:, None] * d[None, :] + 2 * (d * d)[:, None] * d[None, :]


def T_slice(K, k3, Z):
    """transport of the slice-supported tensor with (i,i,k) entries K_ik (i != k) and (i,i,i) entries k3_i to D21."""
    Dl = K.copy(); np.fill_diagonal(Dl, k3 / 3)
    T = Dl @ Z
    return (Z * Z).T @ T + 2 * (Z * T).T @ Z


def main(lb, nh):
    S = bench.load_set('w1024_d16'); W = bench.weights(S, 0); Tm = S['means'][0]
    L, n, _ = W.shape
    tr = []; fc.run(W, trace=tr, slices=2)
    st = tr[lb]; mu, Sm = st['mu'], st['S']
    v = np.diag(Sm); s = np.sqrt(v); al = mu / s
    P = ndtr(al); p = np.exp(-0.5 * al * al) / np.sqrt(2 * np.pi); w2 = p / s
    W1 = W[lb + 1].astype(np.float64)
    mu_ref = np.zeros((L, n))
    for l in range(1, L):
        mu_ref[l] = Tm[l - 1] @ W[l].astype(np.float64)
    acc = {}
    for h in (0, 1):
        rng = np.random.default_rng([777, h]); done = 0
        while done < nh:
            m = 4096
            a = rng.standard_normal((m, n), dtype=np.float32)
            for l in range(lb + 1):
                z = a @ W[l]; a = np.maximum(z, 0)
            y = z - mu_ref[lb].astype(np.float32)
            t = a - Tm[lb].astype(np.float32)
            x = (y * P.astype(np.float32)) @ W[lb + 1]
            zn = a @ W[lb + 1]; yn = zn - mu_ref[lb + 1].astype(np.float32)
            for key, u in (('x', x), ('yn', yn), ('y', y), ('t', t)):
                q = acc.setdefault((h, key), [0, 0, 0, 0])
                q[0] = q[0] + u.sum(0, dtype=np.float64); q[1] = q[1] + (u.T @ u).astype(np.float64)
                q[2] = q[2] + ((u * u).T @ u).astype(np.float64); q[3] = q[3] + (u ** 3).sum(0, dtype=np.float64)
            done += m
    mc = {}
    for h in (0, 1):
        for key in ('x', 'yn', 'y', 't'):
            s1, s2, s21, s3 = [q / nh for q in acc[(h, key)]]
            K = cen21(s21, s2, s1)
            e2 = np.diag(s2); k3 = s3 - 3 * s1 * e2 + 2 * s1 ** 3
            mc[(h, key)] = (K, k3)
    births = []
    for h in (0, 1):
        D = mc[(h, 'yn')][0]; B0 = mc[(h, 'x')][0]
        Ka, k3a = mc[(h, 't')]; Kz, k3z = mc[(h, 'y')]
        Ka = Ka.copy(); np.fill_diagonal(Ka, 0); Kz = Kz.copy(); np.fill_diagonal(Kz, 0)
        sl = T_slice(Ka - (P[:, None] ** 2) * P[None, :] * Kz, k3a - P ** 3 * k3z, W1)
        births.append((D - B0 - sl, D - B0, sl))
    # model: Wick all-distinct transported = T(Wick full) - T(Wick slices)
    Y = (Sm * P[None, :]) @ W1; Z = W1
    wick = ((Y * Y) * w2[:, None]).T @ Z + 2 * ((Y * Z) * w2[:, None]).T @ Y
    Cd = np.diag(Sm)
    Om = (w2[None, :] * (P[:, None] ** 2) * Sm.T ** 2) + 2 * (w2 * P * Cd)[:, None] * P[None, :] * Sm
    np.fill_diagonal(Om, 0)
    om3 = 3 * w2 * P ** 2 * Cd ** 2
    wick_ad = wick - T_slice(Om, om3, Z)
    model_sl = T_slice(st['sl21'] - (P[:, None] ** 2) * P[None, :] * np.where(np.eye(n) > 0, 0, st['D21'] if st['D21'] is not None else 0),
                       st['sl3'] - P ** 3 * (np.diag(st['D21']) if st['D21'] is not None else 0), Z) if 'sl21' in st else None

    def eps(M, A, B):
        return np.sqrt(max(np.sum((M - A) * (M - B)), 0) / np.sum(A * B))
    def nrm(A, B):
        return np.sqrt(max(np.sum(A * B), 0))
    (adA, bA, slA), (adB, bB, slB) = births
    print(f"target {lb + 1}: |D21| {nrm(mc[(0,'yn')][0], mc[(1,'yn')][0]):.3f} |birth| {nrm(bA, bB):.3f} "
          f"|birth_slices| {nrm(slA, slB):.3f} |birth_ad| {nrm(adA, adB):.3f} |model wick_ad| {np.linalg.norm(wick_ad):.3f}")
    print(f"  eps(model wick_ad vs MC birth_ad) = {eps(wick_ad, adA, adB):.3f};  cos = "
          f"{np.sum(wick_ad * (adA + adB)) / np.linalg.norm(wick_ad) / np.linalg.norm(adA + adB):.3f}")
    if model_sl is not None:
        print(f"  eps(model birth slices vs MC) = {eps(model_sl, slA, slB):.3f}")
    # least-squares scale of the Wick term
    c = np.sum(wick_ad * (adA + adB)) / 2 / np.sum(wick_ad ** 2)
    print(f"  best scale on wick_ad: {c:.3f}, eps after scale {eps(c * wick_ad, adA, adB):.3f}")


if __name__ == '__main__':
    main(int(sys.argv[1]), int(sys.argv[2]))


def detail(lb, nh, k4=False):
    import check_birth as cb
    S = bench.load_set('w1024_d16'); W = bench.weights(S, 0); Tm = S['means'][0]
    L, n, _ = W.shape
    k4f = None
    if k4:
        import mc_analyse as ma
        DD = np.load('data/mc1024_mlp0_N262144.npz')
        k4f = {}
        for l in range(1, lb + 1):
            cA = ma.central(DD, 'A', l); cB = ma.central(DD, 'B', l)
            k4f[l] = (0.5 * (cA['K22'] + cB['K22']), 0.5 * (cA['K31'] + cB['K31']), 0.5 * (cA['k4'] + cB['k4']))
    tr = []; fc.run(W, trace=tr, slices=2, k4f=k4f)
    st = tr[lb]; mu, Sm = st['mu'], st['S']
    v = np.diag(Sm); s = np.sqrt(v); P = ndtr(mu / s)
    Z = W[lb + 1].astype(np.float64)
    D = np.load('data/mc1024_mlp0_N65536_a.npz')
    def cen(h, key2, key21, key3, l):
        d = D[f'{h}_s1{key2}'][l]; E2 = D[f'{h}_s2{key2}'][l].astype(float); E21 = D[f'{h}_{key21}'][l].astype(float); E3 = D[f'{h}_{key3}'][l]
        e2 = np.diag(E2); K = cen21(E21, E2, d); k3 = E3 - 3 * d * e2 + 2 * d ** 3
        np.fill_diagonal(K, 0); return K, k3
    def eps(M, A, B):
        return np.sqrt(max(np.sum((M - A) * (M - B)), 0) / max(np.sum(A * B), 1e-30))
    zero = np.zeros(n); Zm = np.zeros((n, n))
    mod21 = st['sl21'].copy(); np.fill_diagonal(mod21, 0)
    Dz = st['D21'].copy(); np.fill_diagonal(Dz, 0)
    out = {}
    for h in 'AB':
        Ka, k3a = cen(h, 't', 's21a', 's3a', lb); Kz, k3z = cen(h, 'y', 's21', 's3', lb)
        out[h] = dict(K=T_slice(Ka, zero, Z), k=T_slice(Zm, k3a, Z), Kz=T_slice((P[:, None] ** 2) * P[None, :] * Kz, zero, Z),
                      kz=T_slice(Zm, P ** 3 * k3z, Z))
    mod = dict(K=T_slice(mod21, zero, Z), k=T_slice(Zm, st['sl3'], Z), Kz=T_slice((P[:, None] ** 2) * P[None, :] * Dz, zero, Z),
               kz=T_slice(Zm, P ** 3 * np.diag(st['D21']), Z))
    for key in mod:
        A, B = out['A'][key], out['B'][key]
        print(key, f"|MC| {np.sqrt(max(np.sum(A*B),0)):.3f} |model| {np.linalg.norm(mod[key]):.3f} eps {eps(mod[key], A, B):.3f}")
    A = out['A']['K'] + out['A']['k'] - out['A']['Kz'] - out['A']['kz']; B = out['B']['K'] + out['B']['k'] - out['B']['Kz'] - out['B']['kz']
    M = mod['K'] + mod['k'] - mod['Kz'] - mod['kz']
    print('birth slices', f"|MC| {np.sqrt(max(np.sum(A*B),0)):.3f} eps {eps(M, A, B):.3f}")
