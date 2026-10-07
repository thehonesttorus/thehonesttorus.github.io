# The output-relevant direction of the D21 error (note XXXVI section 3h).
#   python gradd21.py NET MCPREFIX FREE0 FREE1 ONE0 ONE1
# The output reads the z-level D21 slice at layer k only through the covariance of y_k (facet coefficients) and the
# off-diagonal second moment of z_(k+1): f_k(Gamma) = <prop(k+1, coef_var * diag(W Rk(Gamma) W^T)), e>,
# Rk(Gamma)_ab = (Gamma_ab/2) rho_a Phi_b + (Gamma_ba/2) Phi_a rho_b, whose gradient is
#     grad_k = (rho Phi^T) o (W_(k+1)^T diag(g) W_(k+1)),   g = coef_var(k+1) * (J^T e).
# Evaluated on the free-running D21 errors (FREE0 baseline, FREE1 repaired) and on the one-step errors (own D21 in
# the all-oracle dumps ONE0, ONE1), with e the baseline's output error; % of the baseline n x MSE.
import sys, math, numpy as np
net, pre, f0p, f1p, o0p, o1p = int(sys.argv[1]), sys.argv[2], *sys.argv[3:7]
Wcol = np.load(f"../official/W_off{net}.npy").astype(np.float64); L, n, _ = Wcol.shape
mt = np.load(f"../official/truth_off{net}.npz")["m"].astype(np.float64)
F = np.load(f"{pre}_full.npz")
mu, var = F["mu"].astype(np.float64), F["var"].astype(np.float64)
sig = np.sqrt(var); al = mu / sig; ph = np.exp(-al * al / 2) / math.sqrt(2 * math.pi); Ph = 0.5 * (1 + np.vectorize(math.erf)(al / math.sqrt(2)))
cvar = ph / (2 * sig)
def offd(X):
    X = np.array(X, dtype=np.float64); np.fill_diagonal(X, 0.0); return X
A = {t: np.load(p) for t, p in (("f0", f0p), ("f1", f1p), ("o0", o0p), ("o1", o1p))}
g_ = lambda t, k: A[t][k].astype(np.float64) if k in A[t].files else None
e = g_("f0", f"pk1v_{L - 1}") - mt[-1]; E = float(e @ e)
# adjoint: lam_l = J(L-1 <- l)^T e, lam_(l) = W_(l+1)^T (Ph_(l+1) * lam_(l+1))
lam = {L - 1: e.copy()}
for l in range(L - 2, 0, -1):
    lam[l] = Wcol[l + 1].T @ (Ph[l + 1] * lam[l + 1])
print(f"net {net}: projection of the D21 error on the output-relevant direction, % of baseline n x MSE")
print("    layer k: free base  free rep  (delta) | one-step base  one-step rep  (delta) || L2 err free base rep | one-step base rep | cos(grad, one-step delta)")
tot = np.zeros(6)
for k in range(2, L - 1):
    W = Wcol[k + 1]; gk = cvar[k + 1] * lam[k + 1]
    M = W.T @ (W * gk[:, None])
    r_ = ph[k] / sig[k]
    G = np.outer(r_, Ph[k]) * M; np.fill_diagonal(G, 0.0)
    Dt = offd(F["D21"][k])
    vals = []; l2 = []
    for t, key in (("f0", f"D21_{k}"), ("f1", f"D21_{k}"), ("o0", f"D21own_{k}"), ("o1", f"D21own_{k}")):
        X = g_(t, key)
        if X is None:
            vals.append(np.nan); l2.append(np.nan); continue
        ee = offd(X) - Dt
        vals.append(100 * float(np.sum(G * ee)) / E); l2.append(np.linalg.norm(ee) / np.linalg.norm(Dt))
    dO = offd(g_("o1", f"D21own_{k}") - g_("o0", f"D21own_{k}")) if g_("o0", f"D21own_{k}") is not None and g_("o1", f"D21own_{k}") is not None else None
    cos = float(np.sum(G * dO) / (np.linalg.norm(G) * np.linalg.norm(dO))) if dO is not None else np.nan
    row = [vals[0], vals[1], vals[1] - vals[0], vals[2], vals[3], vals[3] - vals[2]]
    tot += np.nan_to_num(np.array(row))
    print(f"    {k:2d}: {row[0]:+6.2f} {row[1]:+6.2f} ({row[2]:+5.2f}) | {row[3]:+6.2f} {row[4]:+6.2f} ({row[5]:+5.2f}) || "
          f"{l2[0]:.4f} {l2[1]:.4f} | {l2[2]:.4f} {l2[3]:.4f} | {cos:+.3f}", flush=True)
print("    total: " + " ".join(f"{x:+.2f}" for x in tot))
