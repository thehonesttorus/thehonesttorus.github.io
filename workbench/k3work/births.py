# Are the non-gain parts of the fourth-cumulant slices the one-step tree births of the last nonlinearity?
#   python births.py NET   (needs mc2_off{NET}_full.npz, chain_off{NET}_fk_base.npz, ../official/W_off{NET}.npy)
# For jointly Gaussian pre-activations d (covariance C) and y_a = sum_k (c_ak / k!) He_k(d_a), the joint cumulant of
# four linear reads of y is a sum over connected multigraphs on the four reads (no loops), each tree (three covariance
# factors) entering with coefficient 1:
#   star  centred at read t (Hermite order 3 there, 1 at the leaves):  sum_a w_t(a) c_a3 prod_(s != t) M_(a, s)
#   path  s1 - s2 - s3 - s4 (orders 1, 2, 2, 1):  sum_(a, b) w_s2(a) c_a2 M_(a, s1) C_ab w_s3(b) c_b2 M_(b, s4)
# with M_(a, s) = (C (c_1 o w_s))_a. Reads (W_i, W_i, W_i, W_j) give the (3,1) slice K31_ij = kappa(z_i, z_i, z_i, z_j)
# of the next pre-activation z' = W y, and (W_i)^4 its kappa_4 diagonal:
#   K31_tree = 3 S1 + S2 + 6 P1 + 6 P2,   g4_tree = 4 diag(S2) + 12 rowsum(V o X)
#   S1_ij = sum_a W_ia c_a3 M_ai^2 M_aj,  S2_ij = sum_a W_ja c_a3 M_ai^3,  X_ia = W_ia c_a2 M_ai,  V = X C,
#   P1_ij = sum_b V_ib W_ib c_b2 M_bj,    P2_ij = sum_b V_ib M_bi W_jb c_b2,   M = C diag(c1) W^T.
# Evaluated at the true mean and covariance of layer l (Gaussian reference), compared with the true slices of layer
# l + 1 jointly with the gain shapes (3 g var_i C_ij, 3 g var^2), plain and in the read metric, and against the
# chain's own error on the same slices.
import sys, numpy as np
from math import sqrt, pi, erf

net = int(sys.argv[1])
T0 = np.load(f"mc2_off{net}_full.npz"); T = {k: T0[k] for k in T0.files}
ch = np.load(f"chain_off{net}_fk_base.npz")
Wc = np.load(f"../official/W_off{net}.npy").astype(np.float64)          # z_l = Wc[l] @ y_(l-1)
phi = lambda a: np.exp(-0.5 * a * a) / sqrt(2 * pi)
Phi = np.vectorize(lambda a: 0.5 * (1 + erf(a / sqrt(2))))
f64 = lambda a: np.asarray(a, np.float64)


def offd(M):
    M = M.copy(); np.fill_diagonal(M, 0.0); return M


def coefs(l):
    mu, var = f64(T["mu"][l]), f64(T["var"][l]); s = np.sqrt(var); al = mu / s
    return mu, var, Phi(al), phi(al) / s, -al * phi(al) / var, (al * al - 1) * phi(al) / (var * s)


def fit(X, shapes, w):
    # least squares X ~ sum_k b_k S_k in metric w; returns amplitudes and the explained share
    A = np.stack([(S * w).ravel() for S in shapes], 1); y = (X * w).ravel()
    b, *_ = np.linalg.lstsq(A, y, rcond=None)
    return b, 1.0 - float(np.sum((y - A @ b) ** 2)) / float(np.sum(y * y))


for l in range(1, 15):
    mu, var, c1, c2, c3, _ = coefs(l)
    C = f64(T["cov"][l]); W = Wc[l + 1]
    M = C @ (c1[:, None] * W.T)                          # M_ai = (C diag(c1) W^T)_ai = Cov(z_a, u_i)
    Mt = M.T                                             # Mt_ia = M_ai
    S1 = (W * c3[None, :] * Mt * Mt) @ M
    S2 = (Mt ** 3 * c3[None, :]) @ W.T
    X = W * c2[None, :] * Mt; V = X @ C
    P1 = (V * W * c2[None, :]) @ M
    P2 = (V * Mt) @ (W * c2[None, :]).T
    K31b = offd(3 * S1 + S2 + 6 * P1 + 6 * P2)
    g4b = 4 * np.diag(S2) + 12 * np.sum(V * X, axis=1)
    # truth at l + 1 and its gain shapes
    mu1, var1, d1, d2, d3, d4 = coefs(l + 1)
    C1 = f64(T["cov"][l + 1]); K31 = offd(f64(T["K31"][l + 1])); k4 = f64(T["k4"][l + 1])
    Sg31 = offd(3 * var1[:, None] * C1); Sg4 = 3 * var1 * var1
    w31 = d3[:, None] * d1[None, :]; w4 = d4
    one = np.ones_like(K31)
    out = [f"layer {l:2d}->{l + 1:2d}"]
    for nm, Xt, Sg, Bt, w in (("K31", K31, Sg31, K31b, w31), ("g4", k4, Sg4, g4b, w4)):
        for tag, ww in (("plain", np.ones_like(Xt)), ("read", w)):
            (bg,), eg = fit(Xt, [Sg], ww)
            (bb,), eb = fit(Xt, [Bt], ww)
            (bg2, bb2), e2 = fit(Xt, [Sg, Bt], ww)
            R = (Xt - bg * Sg) * ww; Br = Bt * ww
            cr = float(np.sum(R * Br)) / max(np.linalg.norm(R) * np.linalg.norm(Br), 1e-300)
            out.append(f"  {nm} {tag:5s}: gain g {bg:.4f} expl {eg:.3f} | tree b {bb:.3f} expl {eb:.3f} | joint g {bg2:.4f} "
                       f"b {bb2:.3f} expl {e2:.3f} | corr(tree, gain residual) {cr:+.3f}")
    # the chain's own (3,1) slice and kappa_4 diagonal against truth, and the tree's share of the chain's error
    for nm, key, Xt, Bt, w, tr in (("K31", "wk431", K31, K31b, w31, True), ("g4", "g4row", k4, g4b, w4, False)):
        k = f"{key}_{l + 1}"
        if k in ch.files:
            Y = f64(ch[k]); Y = offd(Y.T) if tr else Y
            E = (Xt - Y) * w; Br = Bt * w
            b = float(np.sum(E * Br)) / float(np.sum(Br * Br))
            out.append(f"  chain {nm}: read rel err {np.linalg.norm(E) / np.linalg.norm(Xt * w):.3f}; tree b {b:.3f} explains "
                       f"{1 - float(np.sum((E - b * Br) ** 2)) / float(np.sum(E * E)):.3f} of the (truth - chain) error")
    print("\n".join(out), flush=True)
