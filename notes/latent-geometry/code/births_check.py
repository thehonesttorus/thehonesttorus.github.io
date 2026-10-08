# Check the tree-birth formulas of births.py against Monte Carlo for y_a = c1 d_a + (c2/2) He2(d_a) + (c3/6) He3(d_a)
# (Hermite coefficients known exactly), d ~ N(0, C), z' = W y. Paths: c3 = 0, where kappa_4 = paths + 4-cycles exactly.
# Stars: the part of kappa_4 odd in c3 (c3 -> -c3) is the star to first order in c3.
import numpy as np
rng = np.random.default_rng(1)
n, m = 4, 3
A = rng.standard_normal((n, n)); C = A @ A.T / n + 0.5 * np.eye(n); L = np.linalg.cholesky(C); s2 = np.diag(C)
W = rng.standard_normal((m, n))
def trees(c1, c2, c3):
    M = C @ (c1[:, None] * W.T); Mt = M.T
    S1 = (W * c3[None, :] * Mt * Mt) @ M; S2 = (Mt ** 3 * c3[None, :]) @ W.T
    X = W * c2[None, :] * Mt; V = X @ C
    P1 = (V * W * c2[None, :]) @ M; P2 = (V * Mt) @ (W * c2[None, :]).T
    return 3 * S1 + S2, 6 * P1 + 6 * P2, 4 * np.diag(S2), 12 * np.sum(V * X, axis=1), X, V
def cycles_g4(c2):
    # 4-cycles (orders 2,2,2,2) of the kappa_4 diagonal: 3 labelled cycles, sum_abcd x_a C_ab x_b C_bc x_c C_cd x_d C_da
    out = []
    for i in range(m):
        D = np.diag(W[i] * c2); P = D @ C; out.append(3 * np.trace(P @ P @ P @ P))
    return np.array(out)
def cycles_31(c2):
    # 4-cycles on reads (i,i,i,j): 3 labelled cycles, all with the same value tr(D_i C D_i C D_i C D_j C)
    out = np.zeros((m, m))
    for i in range(m):
        Di = np.diag(W[i] * c2) @ C
        for j in range(m):
            Dj = np.diag(W[j] * c2) @ C; out[i, j] = 3 * np.trace(Di @ Di @ Di @ Dj)
    return out
def mc(c1, c2, c3, N=int(4e7), B=2_000_000, seed=0):
    r = np.random.default_rng(seed); acc = np.zeros(4); S = {k: np.zeros((m, m)) for k in ("m1", "m11", "m21", "m31", "m4")}
    s1 = np.zeros(m); s2_ = np.zeros(m); s3 = np.zeros(m); s4 = np.zeros(m); cnt = 0
    E11 = np.zeros((m, m)); E21 = np.zeros((m, m)); E31 = np.zeros((m, m))
    while cnt < N:
        d = L @ r.standard_normal((n, B)); sd = np.sqrt(s2)[:, None]; x = d / sd
        y = c1[:, None] * d + (c2[:, None] / 2) * sd ** 2 * (x * x - 1) + (c3[:, None] / 6) * sd ** 3 * (x ** 3 - 3 * x)
        u = W @ y; u = u - u.mean(1, keepdims=True) if cnt == 0 else u
        if cnt == 0: m0 = (W @ y).mean(1, keepdims=True)
        u = W @ y - m0
        s1 += u.sum(1); s2_ += (u * u).sum(1); s3 += (u ** 3).sum(1); s4 += (u ** 4).sum(1)
        E11 += u @ u.T; E21 += (u * u) @ u.T; E31 += (u ** 3) @ u.T; cnt += B
    md = s1 / cnt; e2 = s2_ / cnt; e3 = s3 / cnt; e4 = s4 / cnt
    var = e2 - md * md; k4 = e4 - 4 * md * e3 + 6 * md * md * e2 - 3 * md ** 4 - 3 * var * var
    E11 /= cnt; E21 /= cnt; E31 /= cnt
    A_, B_ = md[:, None], md[None, :]; e3i = e3[:, None]; e2i = e2[:, None]
    cov = E11 - A_ * B_
    Ex3y = E31 - B_ * e3i - 3 * A_ * E21 + 3 * A_ * B_ * e2i + 3 * A_ * A_ * E11 - 3 * A_ ** 3 * B_
    K31 = Ex3y - 3 * var[:, None] * cov
    return k4, K31
c1 = rng.uniform(0.3, 1.0, n); c2 = rng.uniform(-0.8, 0.8, n); z = np.zeros(n)
# paths: c3 = 0  ->  exact kappa_4 = trees(paths) + 4-cycles
st31, pa31, stg, pag, X, V = trees(c1, c2, z)
k4mc, K31mc = mc(c1, c2, z)
print("paths, kappa_4 diagonal: MC", np.round(k4mc, 4), " paths + cycles", np.round(pag + cycles_g4(c2), 4))
off = ~np.eye(m, dtype=bool)
print("paths, (3,1) slice offdiag: MC", np.round(K31mc[off], 4), "\n    paths + cycles", np.round((pa31 + cycles_31(c2))[off], 4))
# stars: odd part in c3, with c2 = 0 (no paths), small c3
c3 = rng.uniform(-0.3, 0.3, n); z2 = np.zeros(n)
stp, _, sgp, _, _, _ = trees(c1, z2, c3)
kp, Kp = mc(c1, z2, c3, seed=5); km, Km = mc(c1, z2, -c3, seed=5)
print("stars, kappa_4 diagonal: MC odd part", np.round((kp - km) / 2, 4), " star", np.round(sgp, 4))
print("stars, (3,1) offdiag: MC odd part", np.round(((Kp - Km) / 2)[off], 4), "\n    star", np.round(stp[off], 4))
