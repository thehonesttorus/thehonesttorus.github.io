# Which input drives the off-diagonal second-moment change (checkpoint J steps 3-4, note XXXVI section 3h).
#   python covroute.py NET MCPREFIX DUMP0 DUMP1        (free-running dumps with KEEP_EXTRA=pk1v,K2v,K11)
# For z_l: Delta S2_off(i) = sum_(a != b) W_ia W_ib Delta(C_ab + m_a m_b), y = y_(l-1). Split into
#   mean products   sum_(a != b) W W Delta(m_a m_b)
#   covariance      sum_(a != b) W W Delta C_ab, attributed to the inputs of the pair program at z_(l-1) by the gated
#                   transport coefficients at the true reference (note XXXIV), first order in each input change:
#     C   : Phi_a Phi_b dC_ab + rho_a rho_b C_ab dC_ab                      (Mehler orders 1 and 2)
#     D21 : (dG_ab / 2) rho_a Phi_b + (dG_ba / 2) Phi_a rho_b,             G_ab = kappa(z_a, z_a, z_b)
#     K22 : (dK_ab / 4) rho_a rho_b,                                        K_ab = kappa(z_a, z_a, z_b, z_b)
#     K31 : (dB_ab / 6) c13_a Phi_b + (dB_ba / 6) Phi_a c13_b,             B_ab = kappa(z_a, z_a, z_a, z_b)
#     unary: (rho_a dmu_a + d_var Phi_a dvar_a) Phi_b C_ab + (a <-> b)
#     residual: the chain's actual Delta C minus the five (its closure's departure from these first-order coefficients)
# Each piece is allocated by <prop(l, coef_var * piece), e0 + e1> as in sbudget.py (% of baseline n x MSE).
import sys, math, numpy as np
net, pre, d0p, d1p = int(sys.argv[1]), sys.argv[2], sys.argv[3], sys.argv[4]
Wcol = np.load(f"../official/W_off{net}.npy").astype(np.float64); L, n, _ = Wcol.shape
mt = np.load(f"../official/truth_off{net}.npz")["m"].astype(np.float64)
F = np.load(f"{pre}_full.npz"); A0, A1 = np.load(d0p), np.load(d1p)
mu, var = F["mu"].astype(np.float64), F["var"].astype(np.float64)
sig = np.sqrt(var); al = mu / sig; ph = np.exp(-al * al / 2) / math.sqrt(2 * math.pi); Ph = 0.5 * (1 + np.vectorize(math.erf)(al / math.sqrt(2)))
cvar = ph / (2 * sig)
g = lambda A, k, l: A[f"{k}_{l}"].astype(np.float64) if f"{k}_{l}" in A.files else None
def offd(X):
    X = np.array(X, dtype=np.float64); np.fill_diagonal(X, 0.0); return X
def zmean(A, l):
    m = g(A, "pk1v", l - 1); return None if m is None else Wcol[l] @ m
def prop(l, v):
    for k in range(l + 1, L):
        v = Ph[k] * (Wcol[k] @ v)
    return v
out0, out1 = g(A0, "pk1v", L - 1), g(A1, "pk1v", L - 1)
e0, e1 = out0 - mt[-1], out1 - mt[-1]; s = e0 + e1; E0 = float(e0 @ e0)
print(f"net {net}: {d0p} -> {d1p}: MSE change {100 * float(e1 @ e1 - e0 @ e0) / E0:+.2f}%; off-diagonal second-moment entries, % of baseline n x MSE")
print("    layer: meanprod  cov = [C  D21  K22  K31  unary  resid] || resid/|dC| corr(pred, dC)")
tot = {}
for l in range(2, L):
    k = l - 1
    keys = ("K11", "C_off", "D21", "wk4m", "wk431", "var", "pk1v")
    a0 = {x: g(A0, x, k) for x in keys}; a1 = {x: g(A1, x, k) for x in keys}
    if any(v is None for v in list(a0.values()) + list(a1.values())):
        continue
    W = Wcol[l]
    dCy = offd(a1["K11"] - a0["K11"])
    my0, my1 = a0["pk1v"], a1["pk1v"]
    dmm = offd(np.outer(my1, my1) - np.outer(my0, my0))
    # z_(l-1) gate coefficients at the true reference
    P_, r_ = Ph[k], ph[k] / sig[k]; c13 = -al[k] * ph[k] / var[k]; dPdv = -al[k] * ph[k] / (2 * var[k])
    Cz = offd(0.5 * (F["cov"][k] + F["cov"][k].T))
    dC = offd(a1["C_off"] - a0["C_off"]); dG = offd(a1["D21"] - a0["D21"])
    dK = offd(a1["wk4m"] - a0["wk4m"]); dB = offd((a1["wk431"] - a0["wk431"]).T)        # B_ab = kappa(a,a,a,b)
    dmu = zmean(A1, k) - zmean(A0, k) if k >= 1 else np.zeros(n); dva = a1["var"] - a0["var"]
    R = {}
    R["C"] = np.outer(P_, P_) * dC + np.outer(r_, r_) * Cz * dC
    R["D21"] = 0.5 * dG * np.outer(r_, P_) + 0.5 * dG.T * np.outer(P_, r_)
    R["K22"] = 0.25 * dK * np.outer(r_, r_)
    R["K31"] = dB / 6 * np.outer(c13, P_) + dB.T / 6 * np.outer(P_, c13)
    u = r_ * dmu + dPdv * dva
    R["unary"] = offd(np.outer(u, P_) * Cz + np.outer(P_, u) * Cz)
    pred = sum(R.values()); R["resid"] = dCy - pred
    pieces = {"meanprod": dmm, **R}
    al_ = {}
    for name, M in pieces.items():
        q = np.einsum("ia,ia->i", W @ offd(M), W)              # sum_(a != b) W_ia M_ab W_ib
        al_[name] = 100 * float(prop(l, cvar[l] * q) @ s) / E0
        tot[name] = tot.get(name, 0.0) + al_[name]
    off = ~np.eye(n, dtype=bool)
    print(f"    {l:2d}: {al_['meanprod']:+6.2f}  cov = [{al_['C']:+6.2f} {al_['D21']:+6.2f} {al_['K22']:+6.2f} {al_['K31']:+6.2f} "
          f"{al_['unary']:+6.2f} {al_['resid']:+6.2f}] || {np.linalg.norm(R['resid'][off]) / max(np.linalg.norm(dCy[off]), 1e-300):.3f} "
          f"{np.corrcoef(pred[off], dCy[off])[0, 1]:+.3f}", flush=True)
print("    total: " + " ".join(f"{k} {v:+.2f}" for k, v in tot.items()))
