# The remainder of the kappa3 transport law beyond (legs first order + one facet), against its closed-form next terms.
#   python k3rem.py NET MCPREFIX        (MCPREFIX_full.npz from mcstats post, MCPREFIX_legs.npz from mclegs.py)
# rem = R3 - leg1 - facet (D3 slice), with R3 the measured all-distinct class of y transported (exact, same samples).
# Candidates, each with the coefficient the theorem fixes (second-order cross terms: omega1 omega2 prod_v c_v(1, k1 + k2)):
#   GC1: Gamma x C, C on the doubled index:   (Gamma_uv / 2) c_u(1,3) Phi_v Phi_w C_uw
#   GC2: Gamma x C, C on the single index:    (Gamma_uv / 2) rho_u rho_v Phi_w C_vw
#   DS : C^3 double edge + single edge:       (C_uv^2 / 2) c_u(1,3) rho_v Phi_w C_uw
#   TRI: C^3 triangle:                        rho_a rho_b rho_c C_ab C_bc C_ca      (on a subset of neurons: n^4)
# with Gamma_uv = kappa(z_u, z_u, z_v) (D21), c(1,2) = rho = phi/s, c(1,3) = -alpha phi / s^2, all at the true law of z_l.
import sys, math, numpy as np
net, pre = int(sys.argv[1]), sys.argv[2]
F = np.load(f"{pre}_full.npz"); G = np.load(f"{pre}_legs.npz")
Wcol = np.load(f"../official/W_off{net}.npy").astype(np.float64); L, n, _ = Wcol.shape
mu = F["mu"].astype(np.float64); sd = np.sqrt(F["var"].astype(np.float64)); al = mu / sd
Phi = 0.5 * (1 + np.vectorize(math.erf)(al / math.sqrt(2))); phi = np.exp(-al * al / 2) / math.sqrt(2 * math.pi)
rng = np.random.default_rng(3); sub = np.sort(rng.choice(n, 128, replace=False))
def offd(A):
    A = np.array(A, dtype=np.float64); np.fill_diagonal(A, 0.0); return A
def rel(a, b): return float(np.linalg.norm(a) / np.linalg.norm(b))
def expl(r, x): return 1 - rel(r - x, r) ** 2
def fitn(r, Xs):
    Z = np.stack([x for x in Xs], 1); c, *_ = np.linalg.lstsq(Z, r, rcond=None); return c, 1 - rel(r - Z @ c, r) ** 2
print(f"net {net}: kappa3 law remainder (D3 slice) against closed-form next terms; coefficient 1 and fitted")
for l in range(L - 1):
    W = Wcol[l + 1]; H = W * W; Ph = Phi[l]; rho = phi[l] / sd[l]; c13 = -al[l] * phi[l] / sd[l] ** 2
    X = W * Ph[None, :]; XX = X * X; Wr = W * rho[None, :]
    k3z = F["k3"][l].astype(np.float64); Gm = offd(F["D21"][l]); C = offd(0.5 * (F["cov"][l] + F["cov"][l].T))
    t3 = F["k3"][l + 1].astype(np.float64)
    k3y = F["k3_y"][l].astype(np.float64); Dy = offd(F["D21_y"][l])
    pair = (W * H) @ k3y + 3 * np.einsum("ia,ai->i", H, Dy @ W.T)
    R3 = t3 - pair
    u3 = G["u3"][l].astype(np.float64)
    leg1 = u3 - 3 * np.einsum("ia,ai->i", XX, Gm @ X.T) - (XX * X) @ k3z
    CX = C @ X.T                                                     # (C X^T)_ui = sum_w C_uw X_iw
    fac = 3 * np.einsum("ia,ia->i", Wr, CX.T ** 2 - XX @ (C * C))
    rem = R3 - leg1 - fac
    GX = Gm @ X.T
    gc1 = 3 * np.einsum("iu,ui->i", W * c13[None, :], GX * CX - (Gm * C) @ XX.T)
    Y = Wr * CX.T                                                    # Y_iv = Wr_iv (C X^T)_vi
    gc2 = 3 * (np.einsum("iu,ui->i", Wr, Gm @ Y.T) - np.einsum("iu,ui->i", Wr * X, (Gm * C) @ Wr.T))
    Q = Wr @ (C * C)                                                 # Q_iu = sum_v Wr_iv C_uv^2
    ds = 3 * np.einsum("iu,ui->i", W * c13[None, :], Q.T * CX - (C * C * C) @ (Wr * X).T)
    tri = np.zeros(sub.size)
    for k, i in enumerate(sub):
        v = Wr[i]; M = (C * v[None, :]) @ C
        tri[k] = v @ ((C * M) @ v)
    out = [f"|rem|/|R3| {rel(rem, R3):.3f} |R3|/|t3| {rel(R3, t3):.3f}"]
    for nm, x in (("GC1", gc1), ("GC2", gc2), ("DS", ds)):
        b = float(rem @ x / (x @ x))
        out.append(f"{nm} |.|/|rem| {rel(x, rem):.3f} corr {np.corrcoef(rem, x)[0, 1]:+.2f} coef1 {expl(rem, x):+.3f} fit {b:+.2f}")
    s_ = sub
    out.append(f"TRI(sub) |.|/|rem| {rel(tri, rem[s_]):.3f} corr {np.corrcoef(rem[s_], tri)[0, 1]:+.2f} coef1 {expl(rem[s_], tri):+.3f}")
    allc = gc1 + gc2 + ds
    out.append(f"GC1+GC2+DS coef1 {expl(rem, allc):+.3f}; +TRI(sub) coef1 {expl(rem[s_], allc[s_] + tri):+.3f}")
    c, e = fitn(rem, [gc1, gc2, ds]); out.append(f"joint fit coefs {c[0]:+.2f} {c[1]:+.2f} {c[2]:+.2f} -> {e:.3f}")
    print(f"layer {l:2d}->{l + 1:2d}: " + "\n      ".join(out), flush=True)
