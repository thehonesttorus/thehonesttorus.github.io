# First-entry ledger of the coupled mean / second-moment response (checkpoint J steps 2-3; note XXXVI section 3h).
#   python fentry.py NET MCPREFIX DUMP0 DUMP1 [DUMP2 ...]   (free-running dumps with KEEP_EXTRA=pk1v,K2v,K11)
#   env TRUTHS=full (default) or full,h0,h1 (adds noise-free route Grams from the two Monte Carlo halves)
# Observable state: X_k = (m_k, C^y_k), the post-activation mean and full covariance (diag = var_y), and
# Z_k = (mu_k, C^z_k), the pre-activation mean and covariance (diag = var). Linear maps at the Monte Carlo truth:
#   A_k: Z_k <- X_(k-1)   mu = W m,  C^z = W C^y W^T                                                     (exact)
#   B_k: X_k <- Z_k       dm      = Phi dmu + cv dvar                                                     cv = phi / 2s
#                         dvar_y  = 2 m (1 - Phi) dmu + (Phi - 2 m cv) dvar                              (q = v + m^2)
#                         dC^y_ab = (Phi_a Phi_b + rho_a rho_b C_ab) dC_ab + (u_a Phi_b + Phi_a u_b) C_ab
#                         u = rho dmu + (dPhi / dvar) dvar,  rho = phi / s                     (gated transport, XXXIV)
# For a run, Delta = run - truth (official truth for the means, Monte Carlo for the rest). The first entries are
#   N^z_k = Delta Z_k - A_k Delta X_(k-1),   N^y_k = Delta X_k - B_k Delta Z_k,
# and Delta m_(L-1) = sum_k [prop N^z_k + prop N^y_k] holds EXACTLY: N absorbs every nonlinearity, so the ledger closes
# with no remainder (telescoping). N^y_k is split by route with the gated coefficients:
#   k3, k4 : the diagonal cumulant errors, into m (Gram-Charlier), var_y (through q) and C^y (through Phi);
#   D21    : (dG_ab / 2) rho_a Phi_b + (a <-> b),  K22: (dK_ab / 4) rho_a rho_b,  K31: (dB_ab / 6) c13_a Phi_b + (a <-> b);
#   resid  : the rest of N^y (the chain's response minus the reference linearisation, i.e. closure-derivative defect plus
#            second order);  rep: N^z (the representation step; zero if the chain's C^z is W C^y W^T).
# Forward images delta_X (one linear pass per route) give an exact route decomposition e = sum_X delta_X; the full adjoint
# (a_k, Lambda_k) on X_k, with Lambda_(k-1) = W^T [Lambda_off o K + diag(beta)] W, gives the per-layer entries.
# Reports, per run: route shares <delta_X, e>/|e|^2, the route Gram <delta_X, delta_Y>/|e|^2, per-layer entries;
# per comparison 0 -> r: n Delta MSE = sum_X (own_X + cross_X), own_X = |delta_X^r|^2 - |delta_X^0|^2,
# cross_X = <tau_X, sum_(Y != X) (delta_Y^0 + delta_Y^r)>, tau_X = delta_X^r - delta_X^0, and per layer and route the
# allocation <adjoint(e_0 + e_r), N^r - N^0> (where the adverse change first enters).
import sys, os, math, numpy as np
net, pre, dumps = int(sys.argv[1]), sys.argv[2], sys.argv[3:]
truths = os.environ.get("TRUTHS", "full").split(",")
Wcol = np.load(f"../official/W_off{net}.npy").astype(np.float64); L, n, _ = Wcol.shape
mt = np.load(f"../official/truth_off{net}.npz")["m"].astype(np.float64)
ROUTES = ("k3", "k4", "D21", "K22", "K31", "resid", "rep")
eye = np.eye(n, dtype=bool)


def offd(X):
    X = np.array(X, dtype=np.float64); X[eye] = 0.0; return X


def sym(X):
    X = np.asarray(X, dtype=np.float64); return 0.5 * (X + X.T)


class Truth:
    def __init__(self, tag):
        F = np.load(f"{pre}_{tag}.npz")
        self.mu, self.var = F["mu"].astype(np.float64), F["var"].astype(np.float64)
        self.vy = F["var_y"].astype(np.float64)
        self.k3, self.k4 = F["k3"].astype(np.float64), F["k4"].astype(np.float64)
        self.cov, self.cy = F["cov"], F["cov_y"]
        self.D21, self.K22, self.K31 = F["D21"], F["K22"], F["K31"]

    def Cz(self, k):
        return offd(sym(self.cov[k])) + np.diag(self.var[k])

    def Cy(self, k):
        return offd(sym(self.cy[k])) + np.diag(self.vy[k])


T0 = Truth("full")                                   # the linearisation reference is always the full Monte Carlo
sig = np.sqrt(T0.var); al = T0.mu / sig; ph = np.exp(-al * al / 2) / math.sqrt(2 * math.pi)
Ph = 0.5 * (1 + np.vectorize(math.erf)(al / math.sqrt(2)))
cv = ph / (2 * sig); rho = ph / sig; dPdv = -al * ph / (2 * T0.var)
c3 = -al * ph / (6 * T0.var); c4 = (al * al - 1) * ph / (24 * T0.var * sig)
q3 = ph / (3 * sig); q4 = -al * ph / (12 * T0.var)
f3 = (al * al - 1) * ph / (6 * T0.var * sig); f4 = (3 * al - al ** 3) * ph / (24 * T0.var * T0.var)
c13 = -al * ph / T0.var
Cref = [offd(sym(T0.cov[k])) for k in range(L)]
Kref = [np.outer(Ph[k], Ph[k]) + np.outer(rho[k], rho[k]) * Cref[k] for k in range(L)]


def B(k, dmu, dCz):
    """Reference linearisation X_k <- Z_k; returns (dm, dvar_y, dC^y offdiag)."""
    dvar = np.diag(dCz).copy(); m = mt[k]
    dm = Ph[k] * dmu + cv[k] * dvar
    if k == L - 1:
        return dm, None, None
    dv = 2 * m * (1 - Ph[k]) * dmu + (Ph[k] - 2 * m * cv[k]) * dvar
    u = rho[k] * dmu + dPdv[k] * dvar
    dC = Kref[k] * offd(dCz) + offd(np.outer(u, Ph[k]) * Cref[k] + np.outer(Ph[k], u) * Cref[k])
    return dm, dv, dC


def entries_run(D, T):
    """First entries N[k][route] = (Nm, Nv, NC) for y routes, (None, None, NzC) for rep; checks."""
    g = lambda key, k: D[f"{key}_{k}"].astype(np.float64) if f"{key}_{k}" in D.files else None
    N = []; info = []
    for k in range(L):
        W = Wcol[k]
        # chain state
        if k == 0:
            mp, Cyp = np.zeros(n), np.eye(n); mpt, Cypt = np.zeros(n), np.eye(n)
        else:
            mp = g("pk1v", k - 1); Cyp = offd(sym(g("K11", k - 1))) + np.diag(g("K2v", k - 1))
            mpt, Cypt = mt[k - 1], T.Cy(k - 1)
        Co = g("C_off", k); va = g("var", k)
        WCW = W @ Cyp @ W.T
        Cz = (offd(sym(Co)) if Co is not None else offd(WCW)) + np.diag(va if va is not None else np.diag(WCW))
        dXp_m, dXp_C = mp - mpt, Cyp - Cypt
        dmu = W @ dXp_m                                # chain mu_z = W m_(k-1), truth mu_z = W mt_(k-1)
        dCz = Cz - T.Cz(k)
        NzC = dCz - W @ dXp_C @ W.T
        m = g("pk1v", k)
        dm = m - mt[k]
        if k < L - 1:
            Cy = offd(sym(g("K11", k))) + np.diag(g("K2v", k)); dCy = Cy - T.Cy(k)
            dvy, dCo = np.diag(dCy).copy(), offd(dCy)
        Bm, Bv, BC = B(k, dmu, dCz)
        Nm = dm - Bm
        Nv = dvy - Bv if k < L - 1 else None
        NC = dCo - BC if k < L - 1 else None
        R = {}
        x = g("D3", k)
        if x is not None:
            d = x - T.k3[k]
            R["k3"] = (c3[k] * d, None if k == L - 1 else (q3[k] - 2 * mt[k] * c3[k]) * d,
                       None if k == L - 1 else offd(np.outer(f3[k] * d, Ph[k]) * Cref[k] + np.outer(Ph[k], f3[k] * d) * Cref[k]))
        x = g("g4row", k)
        if x is not None:
            d = x - T.k4[k]
            R["k4"] = (c4[k] * d, None if k == L - 1 else (q4[k] - 2 * mt[k] * c4[k]) * d,
                       None if k == L - 1 else offd(np.outer(f4[k] * d, Ph[k]) * Cref[k] + np.outer(Ph[k], f4[k] * d) * Cref[k]))
        if k < L - 1:
            x = g("D21", k)
            if x is not None:
                dG = offd(x) - offd(T.D21[k])
                R["D21"] = (None, None, 0.5 * dG * np.outer(rho[k], Ph[k]) + 0.5 * dG.T * np.outer(Ph[k], rho[k]))
            x = g("wk4m", k)
            if x is not None:
                dK = offd(x) - offd(T.K22[k])
                R["K22"] = (None, None, 0.25 * dK * np.outer(rho[k], rho[k]))
            x = g("wk431", k)
            if x is not None:
                dB = offd(x.T) - offd(T.K31[k])
                R["K31"] = (None, None, dB / 6 * np.outer(c13[k], Ph[k]) + dB.T / 6 * np.outer(Ph[k], c13[k]))
        zero = lambda v: 0.0 if v is None else v
        rm = Nm - sum(zero(R[r][0]) for r in R)
        rv = None if Nv is None else Nv - sum(zero(R[r][1]) for r in R)
        rC = None if NC is None else NC - sum(zero(R[r][2]) for r in R)
        R["resid"] = (rm, rv, rC)
        R["rep"] = (None, None, NzC)
        N.append(R)
        nrm = lambda X: float(np.linalg.norm(X)) if X is not None else float("nan")
        info.append((k, nrm(NzC) / max(nrm(dCz), 1e-300), nrm(rm) / max(nrm(Nm), 1e-300),
                     nrm(rv) / max(nrm(Nv), 1e-300) if Nv is not None else float("nan"),
                     nrm(rC) / max(nrm(NC), 1e-300) if NC is not None else float("nan"),
                     nrm(NC) / max(nrm(dCo), 1e-300) if NC is not None else float("nan"),
                     sorted(r for r in R if r not in ("resid", "rep"))))
    return N, info


def forward(N):
    """Exact route images delta_X at m_(L-1)."""
    out = {}
    st = {X: (np.zeros(n), np.zeros((n, n))) for X in ROUTES}
    for k in range(L):
        W = Wcol[k]
        for X in ROUTES:
            dm, dC = st[X]
            dmu = W @ dm if k > 0 else np.zeros(n)
            dCz = W @ dC @ W.T if k > 0 else np.zeros((n, n))
            src = N[k].get(X)
            if X == "rep" and src is not None:
                dCz = dCz + src[2]
            bm, bv, bC = B(k, dmu, dCz)
            if src is not None and X != "rep":
                if src[0] is not None: bm = bm + src[0]
                if k < L - 1:
                    if src[1] is not None: bv = bv + src[1]
                    if src[2] is not None: bC = bC + src[2]
            st[X] = (bm, None if k == L - 1 else bC + np.diag(bv))
    return {X: st[X][0] for X in ROUTES}


def adjoint(w):
    """(a_k, Lambda_k) on X_k for the functional <w, m_(L-1)>."""
    a, Lam = w.copy(), np.zeros((n, n)); res = [None] * L
    for k in range(L - 1, -1, -1):
        res[k] = (a, Lam)
        lam = np.diag(Lam).copy(); Lo = offd(Lam)
        h = (Lo * Cref[k]) @ Ph[k]
        m = mt[k]
        alpha = a * Ph[k] + lam * 2 * m * (1 - Ph[k]) + 2 * rho[k] * h
        beta = a * cv[k] + lam * (Ph[k] - 2 * m * cv[k]) + 2 * dPdv[k] * h
        Gam = Lo * Kref[k] + np.diag(beta)
        if k > 0:
            a = Wcol[k].T @ alpha; Lam = Wcol[k].T @ Gam @ Wcol[k]
    return res


def zgam(k, a, Lam):
    lam = np.diag(Lam).copy(); Lo = offd(Lam); h = (Lo * Cref[k]) @ Ph[k]; m = mt[k]
    beta = a * cv[k] + lam * (Ph[k] - 2 * m * cv[k]) + 2 * dPdv[k] * h
    return Lo * Kref[k] + np.diag(beta)


def project(N, adj):
    """per layer, per route entry <adjoint, N>."""
    E = np.zeros((L, len(ROUTES)))
    for k in range(L):
        a, Lam = adj[k]; lam = np.diag(Lam); Lo = offd(Lam)
        for j, X in enumerate(ROUTES):
            src = N[k].get(X)
            if src is None:
                continue
            if X == "rep":
                E[k, j] = float(np.sum(zgam(k, a, Lam) * src[2]))
                continue
            v = 0.0
            if src[0] is not None: v += float(a @ src[0])
            if src[1] is not None: v += float(lam @ src[1])
            if src[2] is not None: v += float(np.sum(Lo * src[2]))
            E[k, j] = v
    return E


runs = [np.load(p) for p in dumps]
outs = [D[f"pk1v_{L - 1}"].astype(np.float64) for D in runs]
es = [o - mt[-1] for o in outs]
E0 = float(es[0] @ es[0])
print(f"net {net}: first-entry ledger; truths {truths}; routes {ROUTES}", flush=True)
# truth self-consistency
tc = [np.linalg.norm(T0.Cz(k) - Wcol[k] @ T0.Cy(k - 1) @ Wcol[k].T) / np.linalg.norm(T0.Cz(k)) for k in range(1, L)]
print("  Monte Carlo consistency |C^z_k - W C^y_(k-1) W^T| / |C^z_k|: " + " ".join(f"{x:.1e}" for x in tc), flush=True)
dirs = {f"e{r}": es[r] for r in range(len(runs))}
dirs.update({f"s{r}": es[0] + es[r] for r in range(1, len(runs))})
ADJ = {name: adjoint(w) for name, w in dirs.items()}
imgs = {}; NN = {}
for t in truths:
    T = T0 if t == "full" else Truth(t)
    for r, D in enumerate(runs):
        N, info = entries_run(D, T)
        img = forward(N); imgs[(t, r)] = img
        if t == truths[0]:
            NN[r] = N
            if True:
                print(f"  run {r} ({dumps[r]}): N diagnostics per layer: |N^z|/|dC^z|, |resid|/|N^y| for (m, var_y, C), |N^y_C|/|dC^y_off|; routes present", flush=True)
                for k, z, rm, rv, rC, nC, pres in info:
                    print(f"    {k:2d}: {z:.3f} | {rm:.3f} {rv:.3f} {rC:.3f} | {nC:.3f} | {','.join(pres)}", flush=True)
        tot = sum(img.values())
        print(f"  [{t}] run {r}: |sum_X delta_X - e| / |e| = {np.linalg.norm(tot - es[r]) / np.linalg.norm(es[r]):.2e}", flush=True)
# ---- per run, own-error route decomposition ----
for r in range(len(runs)):
    e = es[r]; Ee = float(e @ e); img = imgs[(truths[0], r)]
    print(f"\nrun {r} ({dumps[r]}): n MSE {Ee:.4e} (MSE {Ee / n:.4e}); route shares <delta_X, e>/|e|^2 and Gram <delta_X, delta_Y>/|e|^2 [{truths[0]}]")
    print("    route : share   | Gram row (" + " ".join(f"{X:>6s}" for X in ROUTES) + ")")
    for X in ROUTES:
        print(f"    {X:6s}: {float(img[X] @ e) / Ee:+.3f}  | " + " ".join(f"{float(img[X] @ img[Y]) / Ee:+6.3f}" for Y in ROUTES))
    if "h0" in truths and "h1" in truths:
        i0, i1 = imgs[("h0", r)], imgs[("h1", r)]
        print("    noise-free Gram (halves, symmetrised):")
        for X in ROUTES:
            print(f"    {X:6s}: " + " ".join(f"{0.5 * float(i0[X] @ i1[Y] + i1[X] @ i0[Y]) / Ee:+6.3f}" for Y in ROUTES))
    Ep = project(NN[r], ADJ[f"e{r}"])
    print(f"    adjoint check: sum of entries {Ep.sum() / Ee:.6f} (should be 1); per-route forward vs adjoint max diff "
          f"{max(abs(Ep[:, j].sum() - float(img[X] @ e)) for j, X in enumerate(ROUTES)) / Ee:.1e}")
    print("    per layer entries, % of |e|^2: " + " ".join(f"{X:>6s}" for X in ROUTES))
    for k in range(L):
        print(f"      {k:2d}: " + " ".join(f"{100 * Ep[k, j] / Ee:+6.2f}" for j in range(len(ROUTES))), flush=True)
# ---- comparisons ----
for r in range(1, len(runs)):
    s = es[0] + es[r]; dM = float(es[r] @ es[r] - E0)
    i0, i1 = imgs[(truths[0], 0)], imgs[(truths[0], r)]
    print(f"\ncomparison {dumps[0]} -> {dumps[r]}: n Delta MSE = {100 * dM / E0:+.2f}% of baseline")
    print("    route : alloc = own + cross  (% of baseline n MSE)")
    tot = 0.0
    for X in ROUTES:
        tau = i1[X] - i0[X]
        alloc = float(tau @ s); own = float(i1[X] @ i1[X] - i0[X] @ i0[X]); cross = alloc - own
        tot += alloc
        print(f"    {X:6s}: {100 * alloc / E0:+6.2f} = {100 * own / E0:+6.2f} + {100 * cross / E0:+6.2f}")
    print(f"    sum {100 * tot / E0:+.2f} (check vs {100 * dM / E0:+.2f})")
    # pair terms: n Delta MSE = sum_X P_XX + sum_(X < Y) (P_XY + P_YX), P_XY = <tau_X, delta_Y^0 + delta_Y^r>; tau is truth
    # independent, so the Monte Carlo noise enters only through delta^0 + delta^r: noise bar = |S(h0) - S(h1)| / 2
    tau = {X: i1[X] - i0[X] for X in ROUTES}
    Ps = {}
    for t in truths:
        a0, ar = imgs[(t, 0)], imgs[(t, r)]
        Ps[t] = {(X, Y): float(tau[X] @ (a0[Y] + ar[Y])) for X in ROUTES for Y in ROUTES}
    hv = "h0" in truths and "h1" in truths
    S = lambda t, X, Y: Ps[t][(X, Y)] if X == Y else Ps[t][(X, Y)] + Ps[t][(Y, X)]
    print("    own and pair terms (P_XX; P_XY + P_YX), % of baseline" + (", +- half-split noise:" if hv else ":"))
    rows = [(X, Y) for i, X in enumerate(ROUTES) for Y in ROUTES[i:]]
    rows.sort(key=lambda p: -abs(S(truths[0], *p)))
    for X, Y in rows:
        v = 100 * S(truths[0], X, Y) / E0
        if abs(v) < 0.3:
            continue
        nb = f" +- {100 * abs(S('h0', X, Y) - S('h1', X, Y)) / 2 / E0:.2f}" if hv else ""
        print(f"      {X + ('' if X == Y else '-' + Y):10s} {v:+6.2f}{nb}")
    A = ADJ[f"s{r}"]
    P = project(NN[r], A) - project(NN[0], A)
    print("    per layer allocation <adjoint(e0 + er), N^r - N^0>, % of baseline:  " + " ".join(f"{X:>6s}" for X in ROUTES) + "   cum")
    cum = 0.0
    for k in range(L):
        cum += P[k].sum()
        print(f"      {k:2d}: " + " ".join(f"{100 * P[k, j] / E0:+6.2f}" for j in range(len(ROUTES))) + f"   {100 * cum / E0:+6.2f}", flush=True)
