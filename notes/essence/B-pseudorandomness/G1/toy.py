"""G1 toy models (team B).
(A) lazy SRW on the torus Z_m^d: eps-rank (operator norm) of the lambda term W^D - Pi against walk length D, for
    d = 1, 2, 3 at N = m^d ~ 4096..; Dixmier-critical (rank * D ~ const) only at d = 2.
(B) products of a independent Ginibre matrices (free multiplicative): participation ratio n/(a+1) (Fuss-Catalan m2 = a+1)
    and the 90 % / 99 % energy ranks r * a against a; gated version (diag(Phi) with Phi ~ U(0,1) gates) for comparison.
(C) frozen frame through a block for the Ginibre product: freeze the top-r right frame of X = Y G_1..G_a0, carry it
    through G_{a0+1}..G_{2 a0}; relative leg error at every step (theorem: constant in expectation)."""
import numpy as np, json
out = {}
# (A)
eps = 1e-2
for d, m in [(1, 4096), (2, 64), (3, 16)]:
    ks = np.meshgrid(*[np.arange(m)] * d, indexing='ij')
    lam = 1 - sum(np.sin(np.pi * k / m) ** 2 for k in ks) / d
    lam = lam.ravel(); lam = lam[1:] if d else lam
    lam = np.sort(np.abs(lam))[::-1][: m ** d - 1]
    rows = {}
    for D in [1, 2, 4, 8, 16, 32, 64, 128, 256, 512]:
        r = int((lam ** D > eps).sum()); rows[D] = dict(rank=r, rank_x_D=r * D, rank_x_sqrtD=round(r * D ** 0.5, 1), rank_x_D15=round(r * D ** 1.5, 1))
    out[f"torus_d{d}_N{m**d}"] = rows
    print(f"torus d={d} N={m**d}:", {D: (v['rank'], v['rank_x_D']) for D, v in rows.items()})
# (B)
rng = np.random.default_rng(0); n = 1024
for gated in (False, True):
    X = np.eye(n); res = {}
    for a in range(1, 17):
        G = rng.standard_normal((n, n)) / np.sqrt(n)
        if gated:
            ph = rng.uniform(0, 1, n); G = (ph[:, None] * G) * np.sqrt(3.0)   # E phi^2 = 1/3, renormalised
        X = X @ G
        s2 = np.linalg.svd(X, compute_uv=False) ** 2; cs = np.cumsum(s2) / s2.sum()
        pr = s2.sum() ** 2 / (s2 ** 2).sum()
        res[a] = dict(pr_x_a1=round(pr * (a + 1) / n, 3), r90_x_a=round((np.searchsorted(cs, .9) + 1) * a / n, 3), r99_x_a=round((np.searchsorted(cs, .99) + 1) * a / n, 3))
    out["ginibre_gated" if gated else "ginibre"] = res
    print("gated" if gated else "plain", {a: (v['pr_x_a1'], v['r90_x_a'], v['r99_x_a']) for a, v in res.items()})
# (C)
Y = rng.standard_normal((n, n)) / np.sqrt(n)
for a0, r in [(2, 512), (4, 256), (8, 128)]:
    X = Y.copy()
    for _ in range(a0):
        X = X @ (rng.standard_normal((n, n)) / np.sqrt(n))
    _, _, Vt = np.linalg.svd(X, full_matrices=False); Q = Vt[:r].T
    C = X @ Q; Fr = Q.T.copy(); errs = []
    for _ in range(a0):
        G = rng.standard_normal((n, n)) / np.sqrt(n)
        X = X @ G; Fr = Fr @ G
        errs.append(round(float(np.linalg.norm(C @ Fr - X) / np.linalg.norm(X)), 4))
    out[f"frozen_a0{a0}_r{r}"] = errs
    print("frozen a0", a0, "r", r, errs)
json.dump(out, open('toy_results.json', 'w'), indent=1)
