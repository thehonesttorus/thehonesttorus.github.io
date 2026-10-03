# Energy profile of activation outputs F_j = h_{L,j}(X), X ~ N(0, I_d), over exact-mass Gaussian hierarchies:
# level k splits u_k = <q_k, X> (Q orthonormal) into b equal-probability quantile bins (fair b-ary Cantor coding).
# Stratification with one cell per level-K cylinder has variance ratio R = E Var(F|cell) / Var F (gain 1/R);
# any last-point-of-contact walk does no better than 1 - 2 max_k e_k/E (theory), so R bounds what the hierarchy offers.
import numpy as np, sys
from scipy.stats import norm
rng = np.random.default_rng(0)
d = n = 256
def net(L, seed):
    r = np.random.default_rng(seed)
    return [r.standard_normal((n, d if l == 0 else n)).astype(np.float32)*np.sqrt(2/(d if l == 0 else n)) for l in range(L)]
def fwd(Ws, X):
    H = X
    for W in Ws: H = np.maximum(H @ W.T, 0)
    return H
def frames(Ws, J, Xp, Fp):
    out = {}
    out['random'] = np.linalg.qr(rng.standard_normal((d, d)))[0]
    # weights only: linearized network (gates replaced by their mean 1/2)
    Jl = Ws[0].astype(np.float64)
    for W in Ws[1:]: Jl = 0.5*W.astype(np.float64) @ Jl
    # per-output frame: first direction = linearized row of output j, then top right singular vectors of Jl
    U, S, Vt = np.linalg.svd(Jl)
    out['lin-weights'] = [np.linalg.qr(np.column_stack([Jl[j], Vt[:16].T, rng.standard_normal((d, d - 17))]))[0] for j in J]
    # pilot (data-driven upper reference): Stein direction E[X F] then top eigvecs of M = E[(XX^T - I) F]
    pil = []
    for j in J:
        f = Fp[:, j] - Fp[:, j].mean()
        s = Xp.T @ f / len(f)
        M = (Xp*f[:, None]).T @ Xp / len(f)
        ev, V = np.linalg.eigh(M); order = np.argsort(-np.abs(ev))
        pil.append(np.linalg.qr(np.column_stack([s, V[:, order[:16]], rng.standard_normal((d, d - 17))]))[0])
    out['pilot-Stein+Hessian'] = pil
    return out
def ratio(F, U, b, K):
    # cell = tuple of quantile bins of U[:, :K]; returns E Var(F|cell)/Var F, bias-corrected
    edges = norm.ppf(np.arange(1, b)/b)
    cell = np.zeros(len(F), dtype=np.int64)
    for k in range(K): cell = cell*b + np.searchsorted(edges, U[:, k])
    C = b**K
    cnt = np.bincount(cell, minlength=C).astype(float); s1 = np.bincount(cell, F, minlength=C); s2 = np.bincount(cell, F*F, minlength=C)
    ok = cnt > 1
    within = ((s2[ok] - s1[ok]**2/cnt[ok])/(cnt[ok] - 1) * cnt[ok]).sum()/cnt[ok].sum()
    return within/F.var()
for L in [int(a) for a in sys.argv[1:]] or [4, 8, 16]:
    Ws = net(L, 100 + L)
    Xp = rng.standard_normal((20000, d)).astype(np.float32); Fp = fwd(Ws, Xp).astype(np.float64)
    J = [j for j in range(n) if Fp[:, j].var() > 1e-6*Fp.var(0).max()][:8]
    fr = frames(Ws, J, Xp.astype(np.float64), Fp)
    N = 1 << 20
    X = rng.standard_normal((N, d)).astype(np.float32); F = fwd(Ws, X).astype(np.float64)
    # antithetic reference: Var of (F(X)+F(-X))/2 *2 / Var F  (per-evaluation ratio)
    Fm = fwd(Ws, -X[:200000]).astype(np.float64)
    anti = np.mean([((F[:200000, j] + Fm[:, j])/2).var()*2/F[:200000, j].var() for j in J])
    print(f"L={L}: antithetic per-evaluation variance ratio {anti:.3f}", flush=True)
    for name, Q in fr.items():
        for (b, K) in [(2, 6), (2, 12), (4, 3), (4, 6), (8, 4)]:
            Rs = []
            for jj, j in enumerate(J):
                q = Q if name == 'random' else Q[jj]
                U = X.astype(np.float64) @ q[:, :K]
                Rs.append(ratio(F[:, j], U, b, K))
            print(f"   frame={name:20s} b={b} K={K:2d} cells={b**K:5d}: stratification variance ratio R mean {np.mean(Rs):.3f} (min {np.min(Rs):.3f}) -> gain {1/np.mean(Rs):.3f}", flush=True)
