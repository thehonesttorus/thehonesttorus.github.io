# The chain's post-activation covariance error against truth, by input route, and its signed share of the output error.
#   python errroute.py NET MCPREFIX DUMP [DUMP ...]          (free-running dumps with KEEP_EXTRA=pk1v,K2v,K11)
# e_Cy = K11 - Cov_true(y_(l-1)) is split with the gated-transport coefficients at the true reference into the errors of
# the pair program's inputs at z_(l-1) (C_off, D21, K22, K31, mean / variance) plus the value residual (what the closure gets
# wrong from true inputs, measured ~0.4-1.4% in the one-step dump). Each piece is transported to the off-diagonal part of
# S2(z_l) and allocated against the run's own output error e: <prop(l, coef_var * piece), e> / <e, e> (% of the run's MSE).
# If two routes carry large allocations of opposite sign, the run's covariance error is a cancellation between them.
import sys, math, numpy as np
net, pre = int(sys.argv[1]), sys.argv[2]
Wcol = np.load(f"../official/W_off{net}.npy").astype(np.float64); L, n, _ = Wcol.shape
mt = np.load(f"../official/truth_off{net}.npz")["m"].astype(np.float64)
F = np.load(f"{pre}_full.npz")
mu, var = F["mu"].astype(np.float64), F["var"].astype(np.float64)
sig = np.sqrt(var); al = mu / sig; ph = np.exp(-al * al / 2) / math.sqrt(2 * math.pi); Ph = 0.5 * (1 + np.vectorize(math.erf)(al / math.sqrt(2)))
cvar = ph / (2 * sig)
def offd(X):
    X = np.array(X, dtype=np.float64); np.fill_diagonal(X, 0.0); return X
def prop(l, v):
    for k in range(l + 1, L):
        v = Ph[k] * (Wcol[k] @ v)
    return v
for dp in sys.argv[3:]:
    A = np.load(dp); g = lambda k, l: A[f"{k}_{l}"].astype(np.float64) if f"{k}_{l}" in A.files else None
    e = g("pk1v", L - 1) - mt[-1]; E = float(e @ e)
    print(f"net {net} {dp}: MSE {E / n:.4e}; covariance error by input route, % of this run's n x MSE (signed, <tau, e>/<e, e>)")
    print("    layer: C  D21  K22  K31  unary  value || meanprod-err || |e_Cy| rel, value/|e_Cy|")
    tot = {}
    for l in range(2, L):
        k = l - 1
        need = {x: g(x, k) for x in ("K11", "C_off", "D21", "wk4m", "wk431", "var", "pk1v")}
        if any(v is None for v in need.values()):
            continue
        W = Wcol[l]
        Cy_t = offd(0.5 * (F["cov_y"][k] + F["cov_y"][k].T)); eCy = offd(need["K11"]) - Cy_t
        P_, r_ = Ph[k], ph[k] / sig[k]; c13 = -al[k] * ph[k] / var[k]; dPdv = -al[k] * ph[k] / (2 * var[k])
        Cz = offd(0.5 * (F["cov"][k] + F["cov"][k].T))
        eC = offd(need["C_off"]) - Cz; eG = offd(need["D21"]) - offd(F["D21"][k])
        eK = offd(need["wk4m"]) - offd(F["K22"][k]); eB = offd(need["wk431"].T) - offd(F["K31"][k])
        mz = Wcol[k] @ g("pk1v", k - 1) if k >= 1 and g("pk1v", k - 1) is not None else np.zeros(n)
        emu = mz - mu[k]; eva = need["var"] - var[k]
        R = {"C": np.outer(P_, P_) * eC + np.outer(r_, r_) * Cz * eC,
             "D21": 0.5 * eG * np.outer(r_, P_) + 0.5 * eG.T * np.outer(P_, r_),
             "K22": 0.25 * eK * np.outer(r_, r_),
             "K31": eB / 6 * np.outer(c13, P_) + eB.T / 6 * np.outer(P_, c13)}
        u = r_ * emu + dPdv * eva
        R["unary"] = offd(np.outer(u, P_) * Cz + np.outer(P_, u) * Cz)
        R["value"] = eCy - sum(R.values())
        my, myt = need["pk1v"], F["mu_y"][k].astype(np.float64)
        R["meanprod"] = offd(np.outer(my, my) - np.outer(myt, myt))
        al_ = {}
        for name, M in R.items():
            q = np.einsum("ia,ia->i", W @ M, W)
            al_[name] = 100 * float(prop(l, cvar[l] * q) @ e) / E
            tot[name] = tot.get(name, 0.0) + al_[name]
        off = ~np.eye(n, dtype=bool)
        print(f"    {l:2d}: {al_['C']:+6.2f} {al_['D21']:+6.2f} {al_['K22']:+6.2f} {al_['K31']:+6.2f} {al_['unary']:+6.2f} {al_['value']:+6.2f} || "
              f"{al_['meanprod']:+6.2f} || {np.linalg.norm(eCy[off]) / np.linalg.norm(Cy_t[off]):.4f} "
              f"{np.linalg.norm(R['value'][off]) / np.linalg.norm(eCy[off]):.3f}", flush=True)
    print("    total: " + " ".join(f"{k} {v:+.2f}" for k, v in tot.items()), flush=True)
