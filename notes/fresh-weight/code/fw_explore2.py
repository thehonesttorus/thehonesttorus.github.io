# Note XLIV section 7 (exploratory): the top collective mode's eigenvalue error, by source.
#   python -I fw_explore2.py LOCDIR MC2DIR CD2DIR NET LADDER_PY
# u = top eigenvector of the true pre-activation covariance at layer s, r = W_s^T u. Reported relative to lambda_1:
#   total u^T E u; its propagated part u^T P u (Gaussian closure of the chain's state minus that of the true state, s-1)
#   and injected part u^T I u; the injected part by first-order Edgeworth class (chain's own slices at its own state
#   minus the Monte Carlo slices at the true state, s-1), r^T (T_chain - T_true) r; and |u . mu_hat_s|.
import sys, numpy as np
_argv = list(sys.argv)
sys.argv = [_argv[5], ".", ".", "0", "1"]
exec(open(_argv[5]).read().split("cd = np.load")[0])
sys.argv = _argv
loc, mcd, cdd, net = _argv[1], _argv[2], _argv[3], int(_argv[4])
TERMS = ("D21", "k3", "K22", "K31", "k4")


def terms(mu, C, k3, k4, D21, K22, K31):        # as in workbench/k3work/var_attrib.py
    v = np.diag(C).copy(); s = np.sqrt(v); al = mu / s
    R = C / np.outer(s, s); np.fill_diagonal(R, 0.0)
    c = ccoef(al, M + 5)
    sa = s[:, None]; sb = s[None, :]
    _, m1 = relu_var(mu, v)
    out = {}
    T = 0.5 * D21 / sa * series(c, c, R, 2, 1, 0); T = T + T.T; np.fill_diagonal(T, 0.0); out["D21"] = T
    T = (k3[:, None] / 6.0) * sa ** -2 * sb * series(c, c, R, 3, 0, 1); T = T + T.T
    T[np.diag_indices_from(T)] = k3 / 6.0 * 2.0 * c[2] / s - 2.0 * m1 * k3 / 6.0 * c[3] / s ** 2; out["k3"] = T
    T = 0.25 * K22 / (sa * sb) * series(c, c, R, 2, 2, 0); np.fill_diagonal(T, 0.0); out["K22"] = T
    T = (K31 / 6.0) * sa ** -2 * series(c, c, R, 3, 1, 0); T = T + T.T; np.fill_diagonal(T, 0.0); out["K31"] = T
    T = (k4[:, None] / 24.0) * sa ** -3 * sb * series(c, c, R, 4, 0, 1); T = T + T.T
    T[np.diag_indices_from(T)] = k4 / 24.0 * 2.0 * c[3] / s ** 2 - 2.0 * m1 * k4 / 24.0 * c[4] / s ** 3; out["k4"] = T
    return out


W = np.load(f"{loc}/W_off{net}.npy")
c2 = np.load(f"{cdd}/chaindump2_{net}.npz")
F = {h: np.load(f"{mcd}/mc2_off{net}_{h}.npz") for h in ("full", "h0", "h1")}
G = {h: {k: F[h][k] for k in ("mu", "cov", "k3", "k4", "D21", "K22", "K31")} for h in F}
out = c2["out"]
sym = lambda A: 0.5 * (A + A.T)
f64 = lambda a: np.asarray(a, dtype=np.float64)


def chainC(s):
    C = sym(f64(c2[f"C_off_{s}"])); np.fill_diagonal(C, f64(c2[f"var_{s}"])); return C


def zdiag(A):
    A = A.copy(); np.fill_diagonal(A, 0.0); return A


print(f"=== network {net}: top-mode eigenvalue error relative to lambda_1 (mean of halves; half difference in [])",
      flush=True)
for s in range(8, 15):
    Ws = f64(W[s]); t = s - 1
    lam, V = np.linalg.eigh(sym(f64(G["full"]["cov"][s]))); u = V[:, -1]; l1 = lam[-1]
    muh = f64(G["full"]["mu"][s]); align = abs(u @ muh) / np.linalg.norm(muh)
    r = Ws.T @ u
    q = lambda A: float(r @ A @ r)
    mu_c = f64(W[t]) @ f64(out[t - 1]); Cc = chainC(t)
    Tc = terms(mu_c, Cc, f64(c2[f"D3_{t}"]), f64(c2[f"g4row_{t}"]), zdiag(f64(c2[f"D21_{t}"])),
               zdiag(sym(f64(c2[f"wk4m_{t}"]))), zdiag(f64(c2[f"wk431_{t}"]).T))
    Kc = KG(mu_c, Cc)
    tot, prop, inj, cls = [], [], [], {k: [] for k in TERMS}
    for h in ("h0", "h1"):
        mu = f64(G[h]["mu"][t]); C = sym(f64(G[h]["cov"][t]))
        Tt = terms(mu, C, f64(G[h]["k3"][t]), f64(G[h]["k4"][t]), zdiag(f64(G[h]["D21"][t])),
                   zdiag(sym(f64(G[h]["K22"][t]))), zdiag(f64(G[h]["K31"][t])))
        Kt = KG(mu, C)
        Et = chainC(s) - sym(f64(G[h]["cov"][s]))
        tot.append(float(u @ Et @ u)); prop.append(q(Kc - Kt)); inj.append(float(u @ Et @ u) - q(Kc - Kt))
        for k in TERMS:
            cls[k].append(q(Tc[k] - Tt[k]))
    f = lambda a: f"{0.5 * (a[0] + a[1]) / l1:+.2e} [{abs(a[0] - a[1]) / l1:.0e}]"
    print(f"  {s:2d}  |u.mu_hat| {align:.3f}  total {f(tot)}  propagated {f(prop)}  injected {f(inj)}", flush=True)
    print("      injected by class: " + "  ".join(f"{k} {f(cls[k])}" for k in TERMS)
          + f"   sum {f([sum(cls[k][i] for k in TERMS) for i in (0, 1)])}", flush=True)
