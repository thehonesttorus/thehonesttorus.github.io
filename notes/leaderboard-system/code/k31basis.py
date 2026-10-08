# Which basis concentrates the (3,1) slice? (note XXXIX; the open question of note XXXVIII section 9.) The true K31 table
# (wk431 convention, zero diagonal) is projected on the top-K subspaces of candidate objects, two-sided where the object
# defines both sides: captured = |P_U K31 P_V|^2 / |K31|^2 (noise-corrected energy). The oracle is K31's own SVD.
#   python k31basis.py NET MCFULL MCH0 MCH1
import sys, numpy as np

net, ff, f0, f1 = int(sys.argv[1]), sys.argv[2], sys.argv[3], sys.argv[4]
F, H0, H1 = np.load(ff), np.load(f0), np.load(f1)
n = F["mu"].shape[1]
KS = (8, 16, 32, 64)


def zd(x):
    x = np.array(x, dtype=np.float64); np.fill_diagonal(x, 0.0); return x


def flat(D):
    S = (D + D.T) / 2; A = (D - D.T) / 2; off = ~np.eye(n, dtype=bool)
    c = S[off].mean(); al = (S.sum(1) - (n - 1) * c) / (n - 2); be = A.sum(1) / n
    return zd(S - c - al[:, None] - al[None, :] + A - (be[:, None] - be[None, :]))


def top_sym(X, k):
    w, V = np.linalg.eigh((X + X.T) / 2); return V[:, np.argsort(-np.abs(w))[:k]]


print(f"net {net}: captured share of the noise-corrected K31 energy, P_U K31 P_V, at K = " + "/".join(map(str, KS)))
for l in (2, 5, 8, 11, 14, 15):
    K = zd(F["K31"][l].T); dn = zd(H0["K31"][l].T) - zd(H1["K31"][l].T)
    E = np.sum(K ** 2) - np.sum(dn ** 2) / 4
    C = zd(F["cov"][l]) + np.diag(F["var"][l].astype(np.float64)); D = zd(F["D21"][l])
    Uk, sk, Vkt = np.linalg.svd(K); Ud, sd_, Vdt = np.linalg.svd(D); Uf, sf, Vft = np.linalg.svd(flat(D))
    res = {}
    for k in KS:
        Vc = top_sym(C, k)
        cand = {"own SVD (oracle)": (Uk[:, :k], Vkt[:k].T), "cov eigvecs": (Vc, Vc),
                "D21 SVD (U,V)": (Ud[:, :k], Vdt[:k].T), "D21 V both sides": (Vdt[:k].T, Vdt[:k].T),
                "flat D21 SVD": (Uf[:, :k], Vft[:k].T), "cov rows, D21 cols": (Vc, Vdt[:k].T)}
        for name, (U, V) in cand.items():
            P = U.T @ K @ V
            res.setdefault(name, []).append(float(np.sum(P ** 2)) / E)
    print(f"{l:2d} | " + " ; ".join(f"{name} " + "/".join(f"{x:.3f}" for x in v) for name, v in res.items()), flush=True)
