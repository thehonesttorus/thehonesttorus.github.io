# The checkpoint's cell comparison (NCG-20261007-I section 10.2) in the chain's own chart (notes XXXV, XXXVI).
#   python cells.py NET MCPREFIX        (chain_off{NET}_o.npz: one-step dump with every oracle on, inputs MCPREFIX_h0;
#                                         truth MCPREFIX_h1, so input and target sampling noise are independent)
# For each transition l -> l+1 the y-level arrays the chain computed from true inputs (K4v = d, K22 = K) are transported
# to the next pre-activation's three slices (diagonal, (2,2), (3,1)) in several declared tensors, all slices from one tensor:
#   chain : the dumped chain output (diag: seed chart at 2I + K4Q; (2,2): pure trace core at 2I; (3,1): lam core)
#   recon : the same rebuilt here from the arrays (validates the replication of the Ritz step and of lam)
#   c00   : seed chart S = D(d) + B(Khat0) with residual core C(2I, N_E), every slice from that tensor
#   c10   : the same with the residual core transported by the actual metric W W^T
#   c01   : exact transport of the literal retained tensor D(d) + B(K) (= cell 11 in this chart)
#   c6    : c01 plus the transported y (3,1) class B^y of note XXXIV (closed form from the true z_l inputs)
# each with the chain's lam core added (lam read from the dump), and c6 also without it. Reported against Monte Carlo truth:
# relative error and correlation per slice (off-diagonal entries for the matrices), THEORY_3's three diagonal terms, and
# the (2,2) cross term 2 c_ab^T K c_ab on a random subset of pairs (omitted from the full-matrix slices).
import sys, math, numpy as np
net, pre = int(sys.argv[1]), sys.argv[2]
ch = np.load(f"chain_off{net}_o.npz"); T1 = np.load(f"{pre}_h1.npz"); I0 = np.load(f"{pre}_h0.npz")
Wcol = np.load(f"../official/W_off{net}.npy").astype(np.float64)
n = Wcol.shape[1]; off = ~np.eye(n, dtype=bool); one = np.ones(n)
cA, cI = 6.0 / (n + 4.0), -3.0 / ((n + 2.0) * (n + 4.0))
P = lambda v: cA * v + cI * v.sum() * one
Phi = lambda x: 0.5 * (1 + np.vectorize(math.erf)(x / math.sqrt(2))); phi = lambda x: np.exp(-x * x / 2) / math.sqrt(2 * math.pi)
rng = np.random.default_rng(0)


def rel(p, t):
    return np.linalg.norm(p - t) / np.linalg.norm(t)


def core(M, N):
    """slices of the normalised product core C(M, N): diag, (2,2), (3,1) with (3,1)_ij = kappa(i,i,i,j)."""
    m, nd = np.diag(M), np.diag(N)
    return m * nd, (m[:, None] * nd[None, :] + m[None, :] * nd[:, None] + 4 * M * N) / 6, (m[:, None] * N + M * nd[:, None]) / 2


def pair_lift(W, H, K, d):
    """slices of W#[D(d) + B(K)] without the (2,2) cross term 2 c_ab^T K c_ab."""
    HK = H @ K
    W3 = W * H
    diag = (H * H) @ d + 3 * np.einsum("ia,ia->i", HK, H)
    k22 = HK @ H.T + (H * d[None, :]) @ H.T
    k31 = 3 * (HK * W) @ W.T + (W3 * d[None, :]) @ W.T
    return diag, k22, k31


def s31_lift(W, H, B):
    G = B @ W.T; X = W * G.T; W3 = W * H
    return 4 * np.einsum("ia,ai->i", W3, G), 2 * (H @ X.T + (H @ X.T).T), W3 @ G + 3 * (H * G.T) @ W.T


def gate(mu, var):
    s = np.sqrt(var); a = mu / s; Ph = Phi(a); ph = phi(a)
    m = mu * Ph + s * ph; v = (mu * mu + var) * Ph + mu * s * ph - m * m
    dl, dl1 = ph / s, -a * ph / var
    c33 = 6 * Ph - 6 * m * dl + 3 * m * m * dl1
    return dict(Ph=Ph, dl=dl, dl1=dl1, e31=3 * (v - m * m) * (1 - Ph), e32=6 * m * (1 - Ph) + 3 * (m * m - v) * dl,
                e33=c33 - 3 * v * dl1)


print(f"net {net}: inputs {pre}_h0 (true z_l), truth {pre}_h1 (z_(l+1)); rel err / corr per slice (diag | (2,2) | (3,1))")
for l in range(1, 15):
    need = [f"K22_{l}", f"K4v_{l}", f"K2v_{l}", f"g4rowown_{l + 1}", f"wk4mown_{l + 1}", f"wk431own_{l + 1}",
            f"C_off_{l + 1}", f"var_{l + 1}", f"C_off_{l}", f"var_{l}", f"D21_{l}", f"wk4m_{l}", f"wk431_{l}"]
    if any(k not in ch.files for k in need):
        continue
    W = Wcol[l + 1]; H = W * W
    K = ch[f"K22_{l}"].astype(np.float64); K = 0.5 * (K + K.T); np.fill_diagonal(K, 0.0)
    d = ch[f"K4v_{l}"].astype(np.float64); vy = ch[f"K2v_{l}"].astype(np.float64)
    # Ritz approximation exactly as the code: sketch K W[:p]^T, one power iteration, p = K4Q_RANK + 8 = 12 pairs
    p_ = 12
    Q, _ = np.linalg.qr(K @ W[:p_].T); Q, _ = np.linalg.qr(K @ (K @ Q))
    lb, Ub = np.linalg.eigh(Q.T @ K @ Q); V = Q @ Ub
    Kh0 = (V * lb) @ V.T; np.fill_diagonal(Kh0, 0.0)
    E = K - Kh0; r = E @ one; pv = P(r); u = d + K @ one
    M = W @ W.T; NE = (W * pv[None, :]) @ W.T; Nfull = (W * P(u)[None, :]) @ W.T
    # lam core, from the dump: wk431own = lam C_off(z_(l+1)) with the oracle covariance; s_off^2 with the oracle variance
    Cz = ch[f"C_off_{l + 1}"].astype(np.float64); vz = ch[f"var_{l + 1}"].astype(np.float64)
    own431 = ch[f"wk431own_{l + 1}"].astype(np.float64)
    lam = float(np.sum(own431[off] * Cz[off]) / np.sum(Cz[off] ** 2))
    so2 = vz - H @ vy; Czf = Cz.copy(); np.fill_diagonal(Czf, 0.0)
    lam_core = (2 * lam * so2, lam * (so2[:, None] + so2[None, :]) / 3, lam * Czf)
    # cells (pair part)
    LS = pair_lift(W, H, Kh0, d); LE = pair_lift(W, H, E, np.zeros(n))
    c2I, cM = core(2 * np.eye(n), NE), core(M, NE)
    c00 = tuple(a + b for a, b in zip(LS, c2I)); c10 = tuple(a + b for a, b in zip(LS, cM))
    c01 = tuple(a + b for a, b in zip(LS, LE))
    chain_pair = (c00[0], (np.diag(Nfull)[:, None] + np.diag(Nfull)[None, :]) / 3, np.zeros((n, n)))
    # the y (3,1) class from the true z_l inputs (note XXXIV), transported
    g = gate(I0["mu"][l].astype(np.float64), ch[f"var_{l}"].astype(np.float64))
    C0 = ch[f"C_off_{l}"].astype(np.float64).copy(); np.fill_diagonal(C0, 0.0)
    Gm = ch[f"D21_{l}"].astype(np.float64).copy(); np.fill_diagonal(Gm, 0.0)
    Kz = ch[f"wk4m_{l}"].astype(np.float64).copy(); np.fill_diagonal(Kz, 0.0)
    Bz = ch[f"wk431_{l}"].astype(np.float64).T.copy(); np.fill_diagonal(Bz, 0.0)     # Bz_ab = kappa(a,a,a,b)
    Ph, dl, dl1 = g["Ph"], g["dl"], g["dl1"]
    By = (g["e31"][:, None] * Ph[None, :] * C0 + 0.5 * Gm * g["e32"][:, None] * Ph[None, :]
          + 0.5 * Gm.T * g["e31"][:, None] * dl[None, :] + 0.25 * Kz * g["e32"][:, None] * dl[None, :]
          + Bz / 6 * g["e33"][:, None] * Ph[None, :] + Bz.T / 6 * g["e31"][:, None] * dl1[None, :])
    np.fill_diagonal(By, 0.0)
    L31 = s31_lift(W, H, By)
    c6 = tuple(a + b for a, b in zip(c01, L31))
    # truth
    t = (T1["k4"][l + 1].astype(np.float64), T1["K22"][l + 1].astype(np.float64), T1["K31"][l + 1].astype(np.float64))
    nz = (np.linalg.norm(I0["k4"][l + 1] - T1["k4"][l + 1]) / math.sqrt(2) / np.linalg.norm(t[0]),
          np.linalg.norm((I0["K22"][l + 1] - T1["K22"][l + 1])[off]) / math.sqrt(2) / np.linalg.norm(t[1][off]),
          np.linalg.norm((I0["K31"][l + 1] - T1["K31"][l + 1])[off]) / math.sqrt(2) / np.linalg.norm(t[2][off]))
    own = (ch[f"g4rowown_{l + 1}"].astype(np.float64), ch[f"wk4mown_{l + 1}"].astype(np.float64),
           ch[f"wk431own_{l + 1}"].astype(np.float64).T)
    plus = lambda X: tuple(a + b for a, b in zip(X, lam_core))

    def row(X):
        out = []
        for k, (p, tt) in enumerate(zip(X, t)):
            pp, tv = (p, tt) if k == 0 else (p[off], tt[off])
            out.append(f"{rel(pp, tv):.3f}/{np.corrcoef(pp, tv)[0, 1]:+.3f}")
        return " | ".join(out)
    rec = plus(chain_pair)
    vrec = (rel(rec[0], own[0]), rel(rec[1][off], own[1][off]), rel(rec[2][off], own[2][off]))
    # THEORY_3's diagonal ledger of R_W(E) (the chain's diagonal is c01's minus R_W(E), at fixed lam)
    s_ = r.sum(); gg = (r - s_ / (2 * (n - 1)) * one) / (n - 2)
    Ag = np.outer(gg, one) + np.outer(one, gg) - 2 * np.diag(gg); E0 = E - Ag; rH = H @ one
    Dm = (rH - 2) * (H @ pv); Dk = 3 * np.einsum("ia,ia->i", H @ Ag, H) - rH * (H @ pv); Du = 3 * np.einsum("ia,ia->i", H @ E0, H)
    res = t[0] - own[0]
    led = " ".join(f"{nm} {np.linalg.norm(x) / np.linalg.norm(t[0]):.3f} (corr w/ residual {np.corrcoef(x, res)[0, 1]:+.2f})"
                   for nm, x in (("D_metric", Dm), ("D_known", Dk), ("D_unres", Du)))
    # (2,2) cross term on 3000 random pairs
    ia = rng.integers(0, n, 3000); ib = rng.integers(0, n, 3000); keep = ia != ib; ia, ib = ia[keep], ib[keep]
    Cab = W[ia] * W[ib]; cross = 2 * np.einsum("pa,pa->p", Cab @ K, Cab)
    cross_rel = np.linalg.norm(cross) / np.linalg.norm(c01[1][ia, ib])
    print(f"layer {l:2d}->{l + 1:2d}: lam {lam:.5f}; MC noise of truth {nz[0]:.3f} {nz[1]:.3f} {nz[2]:.3f}; replication "
          f"{vrec[0]:.1e} {vrec[1]:.1e} {vrec[2]:.1e}\n"
          f"    chain {row(own)}\n    c00   {row(plus(c00))}\n    c10   {row(plus(c10))}\n"
          f"    c01   {row(plus(c01))}\n    c6    {row(plus(c6))}\n    c6-noλ {row(c6)}\n"
          f"    diag ledger: {led}\n    (3,1)_y transport size vs truth: {np.linalg.norm(L31[0]) / np.linalg.norm(t[0]):.3f} "
          f"{np.linalg.norm(L31[1][off]) / np.linalg.norm(t[1][off]):.3f} {np.linalg.norm(L31[2][off]) / np.linalg.norm(t[2][off]):.3f}; "
          f"(2,2) cross term / main term on 3000 pairs: {cross_rel:.4f}", flush=True)
