# What is in the true higher-cumulant slices beyond the gain, and is the remainder collective?
#   python latent_slices.py NET   (needs mc2_off{NET}_full.npz and chain_off{NET}_fk_base.npz in the working dir)
# Per layer l of Monte Carlo truth (pre-activation z_l):
#   (1) the covariance spectrum: leading eigenvalues over mean variance, trace share of the top K, and the overlap of
#       the mean direction with the leading eigenvectors;
#   (2) the scale-mixture (gain) share of each slice: z = G z~ with Var G^2 = g gives, to first order in g,
#       k3 = 1.5 g mu var, D21_ij = (g/2)(2 mu_i C_ij + mu_j var_i), k4 = 3 g var^2, K22_ij = g (var_i var_j + 2 C_ij^2),
#       K31_ij = 3 g var_i C_ij; g is fitted per slice and the explained share 1 - |X - g S|^2 / |X|^2 is reported,
#       plain and in the read metric the pair program uses (row and column weights by the Hermite coefficients);
#   (3) the residual X - g S: the share of its read energy in the top-K eigenspace of C on the column index;
#   (4) the chain's slices against truth: relative read error, its gain-shaped share and its collective share.
import sys, numpy as np
from math import sqrt, pi, erf

net = int(sys.argv[1])
T0 = np.load(f"mc2_off{net}_full.npz"); ch = np.load(f"chain_off{net}_fk_base.npz")
T = {k: T0[k] for k in T0.files}   # read each array once
phi = lambda a: np.exp(-0.5 * a * a) / sqrt(2 * pi)
Phi = np.vectorize(lambda a: 0.5 * (1 + erf(a / sqrt(2))))
f64 = lambda a: np.asarray(a, np.float64)
print(f"net {net}: truth n = {int(T['n'])}; chain keys e.g. {sorted(ch.files)[:6]}")


def offd(M):
    M = M.copy(); np.fill_diagonal(M, 0.0); return M


def share(X, S, w=None):
    # least-squares amplitude of shape S in X (optionally weighted), and the explained share
    if w is not None:
        X, S = X * w, S * w
    ss = float(np.sum(S * S)); g = float(np.sum(X * S)) / ss if ss > 0 else 0.0
    return g, 1.0 - float(np.sum((X - g * S) ** 2)) / max(float(np.sum(X * X)), 1e-300)


def colshare(X, U):
    # share of X's energy captured by projecting its column index on span(U)
    return float(np.sum((X @ U) ** 2)) / max(float(np.sum(X * X)), 1e-300)


for l in range(1, 16):
    mu, var = f64(T["mu"][l]), f64(T["var"][l]); C = f64(T["cov"][l])
    k3, k4 = f64(T["k3"][l]), f64(T["k4"][l])
    D21, K22, K31 = offd(f64(T["D21"][l])), offd(f64(T["K22"][l])), offd(f64(T["K31"][l]))
    s = np.sqrt(var); al = mu / s
    c1, c2, c3, c4 = Phi(al), phi(al) / s, -al * phi(al) / var, (al * al - 1) * phi(al) / (var * s)
    ev, U = np.linalg.eigh(C); ev, U = ev[::-1], U[:, ::-1]
    mv = var.mean(); muh = mu / np.linalg.norm(mu)
    ov = (U.T @ muh) ** 2
    print(f"\nlayer {l:2d}  alpha sd {al.std():.2f}  eig/mean var: " + " ".join(f"{e / mv:.1f}" for e in ev[:6])
          + "  trace share top " + " ".join(f"{k}:{ev[:k].sum() / ev.sum():.3f}" for k in (1, 4, 16, 64, 256))
          + f"  |<u1,mu^>|^2 {ov[0]:.3f}  top16 {ov[:16].sum():.3f}  top64 {ov[:64].sum():.3f}")
    shapes = {
        "k3": (k3, 1.5 * mu * var, c3),
        "D21": (D21, offd(0.5 * (2 * mu[:, None] * C + var[:, None] * mu[None, :])), (c2[:, None], c1[None, :])),
        "k4": (k4, 3 * var * var, c4),
        "K22": (K22, offd(var[:, None] * var[None, :] + 2 * C * C), (c2[:, None], c2[None, :])),
        "K31": (K31, offd(3 * var[:, None] * C), (c3[:, None], c1[None, :])),
    }
    for nm, (X, S, w) in shapes.items():
        wm = w if not isinstance(w, tuple) else w[0] * w[1]
        g0, e0 = share(X, S); g1, e1 = share(X, S, wm)
        line = f"  {nm:4s} gain fit g {g0:+.4f} explains {e0:+.3f}; read-weighted g {g1:+.4f} explains {e1:+.3f}"
        if X.ndim == 2:
            R = (X - g1 * S) * wm; Xr = X * wm
            line += ("  | read energy in top-K col span: truth " + " ".join(f"{k}:{colshare(Xr, U[:, :k]):.3f}" for k in (16, 64, 256))
                     + "; gain residual " + " ".join(f"{k}:{colshare(R, U[:, :k]):.3f}" for k in (16, 64, 256)))
        print(line, flush=True)
    # the chain's slices against truth (read metric)
    cmp = {"D3": ("D3", k3, 1.5 * mu * var, c3), "g4": ("g4row", k4, 3 * var * var, c4)}
    for nm, (key, X, S, w) in cmp.items():
        k = f"{key}_{l}"
        if k in ch.files:
            Y = f64(ch[k]); E = (Y - X) * w; Xr = X * w; Sr = S * w
            ge = float(np.sum(E * Sr)) / float(np.sum(Sr * Sr))
            print(f"  chain {nm}: read rel err {np.linalg.norm(E) / np.linalg.norm(Xr):.3f}; gain-shaped share of the error "
                  f"{ge * ge * float(np.sum(Sr * Sr)) / float(np.sum(E * E)):.3f}", flush=True)
    for nm, key, X, S, (wr, wc), trans in (("D21", "D21", D21, shapes["D21"][1], shapes["D21"][2], False),
                                          ("K22", "wk4m", K22, shapes["K22"][1], shapes["K22"][2], False),
                                          ("K31", "wk431", K31, shapes["K31"][1], shapes["K31"][2], True)):
        k = f"{key}_{l}"
        if k in ch.files:
            Y = offd(f64(ch[k])); Y = Y.T if trans else Y
            wm = wr * wc; E = (Y - X) * wm; Xr = X * wm; Sr = S * wm
            ge = float(np.sum(E * Sr)) / float(np.sum(Sr * Sr))
            print(f"  chain {nm}: read rel err {np.linalg.norm(E) / np.linalg.norm(Xr):.3f}; gain-shaped share "
                  f"{ge * ge * float(np.sum(Sr * Sr)) / float(np.sum(E * E)):.3f}; error energy in top-K col span "
                  + " ".join(f"{kk}:{colshare(E, U[:, :kk]):.3f}" for kk in (16, 64, 256)), flush=True)
