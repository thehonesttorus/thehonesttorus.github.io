# Continuation-equivalence ("quotient before truncation") tests on the paired arrays (note XXXVI section 3j).
#   python nullalign.py NET MCPREFIX DUMP0 [DUMP1 ...]      (free-running dumps: var, C_off, D3, D21, g4row, wk4m, wk431, pk1v)
# Exact facts used (THEORY.md; verified numerically in nullsrc/verify_null.py):
#  * For every degree-1 homogeneous continuation at a Gaussian reference N(mu, Sigma) the source tuple
#        (dSigma, dk3, dk4) = (G/12, mu.G/4, Sigma.G)     (normalised symmetric products)
#    is invisible; for degree-2 observables (the post-activation covariance) the invisible tuple is (0, mu.G/4, Sigma.G).
#  * The mean-one scale ("gain") mode z = S Z0, E S = 1, Var S = t, has the law tangent
#        (t (Sigma + mu mu^T), 6 t mu.Sigma, 12 t Sigma.Sigma),
#    invisible to every future mean; per unit cv (s^2 + mu^2) + 6 c3 mu s^2 + 12 c4 s^4 = 0 identically.
# Slice templates (B_ab = kappa(a,a,a,b), D21_ab = kappa(a,a,b), K22_ab = kappa(a,a,b,b)):
#    Sigma.G : diag var_a G_aa;  K22 (var_a G_bb + var_b G_aa + 4 C_ab G_ab)/6;  B (var_a G_ab + C_ab G_aa)/2
#    gain k4 = 12 Sigma.Sigma : diag 12 var^2;  K22 4 (var_a var_b + 2 C_ab^2);  B 12 var_a C_ab
#    gain k3 = 6 mu.Sigma : diag 6 mu var;  D21 2 (2 mu_a C_ab + mu_b var_a)
# Reports, per layer, for the truth and each chain run:
#  (1) which metric makes the kappa4 slices a trace core: covariance metric Sigma.G (per-pair least squares) against the
#      Euclidean 2I.G the chain uses (METRIC_C = 2): residual fractions of the (2,2) and (3,1) slices;
#  (2) gain-mode amplitudes t fitted separately on each slice family (kappa4: diag, K22, B; kappa3: diag, D21) and the
#      slice energy they explain - consistency across families is the test of the gain-mode structure;
#  (3) for the chains: the per-unit mean error a gain-amplitude inconsistency predicts, 6 c3 mu var dt3 + 12 c4 var^2 dt4,
#      against the actual query error P3 + P4 = c3 dk3 + c4 dk4 (cosine, explained fraction);
#  (4) alignment of each chain's lower-order errors with the null companions (dG/12, mu.dG/4) of its kappa4 error trace
#      core dG (covariance metric at the true Sigma); for runs >= 1 also of the change against the baseline.
import sys, math, numpy as np
net, pre, dumps = int(sys.argv[1]), sys.argv[2], sys.argv[3:]
Wcol = np.load(f"../official/W_off{net}.npy").astype(np.float64); L, n, _ = Wcol.shape
mt = np.load(f"../official/truth_off{net}.npz")["m"].astype(np.float64)
eye = np.eye(n, dtype=bool); erf = np.vectorize(math.erf)
METRIC_C = 2.0


def offd(X):
    X = np.array(X, dtype=np.float64); X[eye] = 0.0; return X


def sym(X):
    X = np.asarray(X, dtype=np.float64); return 0.5 * (X + X.T)


F = np.load(f"{pre}_full.npz")
TV = {k: F[k] for k in ("mu", "var", "k3", "k4")}
TM = {k: F[k] for k in ("cov", "D21", "K22", "K31")}
runs = [np.load(p) for p in dumps]
g = lambda D, key, k: D[f"{key}_{k}"].astype(np.float64) if f"{key}_{k}" in D.files else None


def fit_cov(dg, K22, B, var, C):
    """Least-squares Sigma.G trace core from (diag, K22, B) slices; G_aa exact from the diagonal, G_ab per pair."""
    Gd = dg / var
    va, vb = var[:, None], var[None, :]
    y1 = K22 - (va * Gd[None, :] + vb * Gd[:, None]) / 6.0          # = (2/3) C_ab x
    y2 = B - C * Gd[:, None] / 2.0                                   # = var_a x / 2
    y3 = B.T - C * Gd[None, :] / 2.0                                 # = var_b x / 2
    a1, a2, a3 = (2.0 / 3.0) * C, va / 2.0, vb / 2.0
    x = (a1 * y1 + a2 * y2 + a3 * y3) / (a1 * a1 + a2 * a2 + a3 * a3)
    x = offd(x)
    G = x + np.diag(Gd)
    r22 = offd(y1 - a1 * x); rB = offd(y2 - a2 * x)
    return G, r22, rB


def fit_euc(dg, K22, B):
    """Euclidean (2I).G core: diag 2 G_aa, K22 (2 G_aa + 2 G_bb)/6, B 2 G_ab / 2."""
    Gd = dg / METRIC_C
    x = offd(B + B.T) / METRIC_C
    r22 = offd(K22 - METRIC_C * (Gd[:, None] + Gd[None, :]) / 6.0)
    rB = offd(B - METRIC_C * x / 2.0)
    return r22, rB


def tfit(slice_, tmpl, mask=None):
    s, t_ = (slice_, tmpl) if mask is None else (slice_[mask], tmpl[mask])
    tt = float(np.sum(t_ * t_)); a = float(np.sum(t_ * s)) / tt if tt > 0 else float("nan")
    ex = 1.0 - float(np.sum((s - a * t_) ** 2)) / max(float(np.sum(s * s)), 1e-300)
    return a, ex


def coef(k):
    mu, var = TV["mu"][k].astype(np.float64), TV["var"][k].astype(np.float64)
    s = np.sqrt(var); a = mu / s; ph = np.exp(-a * a / 2) / math.sqrt(2 * math.pi); Ph = 0.5 * (1 + erf(a / math.sqrt(2)))
    return dict(mu=mu, var=var, s=s, a=a, ph=ph, Ph=Ph, cv=ph / (2 * s), c3=-a * ph / (6 * var), c4=(a * a - 1) * ph / (24 * var * s), rho=ph / s)


rf = lambda X, Y: math.sqrt(float(np.sum(X * X)) / max(float(np.sum(Y * Y)), 1e-300))
print(f"net {net}: continuation-equivalence tests; runs {dumps}", flush=True)
print("  (1)-(2) per layer: kappa4 trace-core residual fractions |r|/|slice| for (2,2) and (3,1), covariance metric | Euclidean 2I;"
      "\n          gain amplitudes t (and explained fraction) from kappa4 diag / K22 / B and kappa3 diag / D21")
for k in range(1, L):
    cf = coef(k); var, mu = cf["var"], cf["mu"]
    C = offd(sym(TM["cov"][k])); Sg = C + np.diag(var)
    tv4 = (12 * var * var, 4 * (np.outer(var, var) + 2 * C * C), 12 * var[:, None] * C)
    tv3 = (6 * mu * var, 2 * (2 * mu[:, None] * C + np.outer(var, mu)))
    rows = []
    srcs = [("truth", TV["k4"][k].astype(np.float64), offd(TM["K22"][k]), offd(TM["K31"][k]), TV["k3"][k].astype(np.float64), offd(TM["D21"][k]))]
    for r, D in enumerate(runs):
        if g(D, "g4row", k) is None or g(D, "wk431", k) is None:
            continue
        srcs.append((f"run{r}", g(D, "g4row", k), offd(g(D, "wk4m", k)), offd(g(D, "wk431", k).T), g(D, "D3", k), offd(g(D, "D21", k))))
    for name, dg, K22, B, d3, D21 in srcs:
        G, r22, rB = fit_cov(dg, K22, B, var, C)
        e22, eB = fit_euc(dg, K22, B)
        t4d, x4d = tfit(dg, tv4[0]); t4p, x4p = tfit(K22, tv4[1], ~eye); t4b, x4b = tfit(B, tv4[2], ~eye)
        t3d, x3d = tfit(d3, tv3[0]); t3p, x3p = tfit(D21, tv3[1], ~eye)
        rows.append(f"      {name:5s}: cov {rf(r22, K22):.3f} {rf(rB, B):.3f} | euc {rf(e22, K22):.3f} {rf(eB, B):.3f} || "
                    f"t4 {t4d:+.2e}({x4d:+.2f}) {t4p:+.2e}({x4p:+.2f}) {t4b:+.2e}({x4b:+.2f}) | t3 {t3d:+.2e}({x3d:+.2f}) {t3p:+.2e}({x3p:+.2f})")
    print(f"    layer {k:2d}:\n" + "\n".join(rows), flush=True)

print("  (3) per layer: the mean error a gain-amplitude inconsistency predicts (6 c3 mu var dt3 + 12 c4 var^2 dt4, dt from the"
      "\n      diagonal fits) against the chain's query error P3 + P4: cos, explained fraction 1 - |P - pred|^2/|P|^2, |pred|/|P|")
for r, D in enumerate(runs):
    line = []
    for k in range(1, L):
        if g(D, "g4row", k) is None:
            continue
        cf = coef(k); var, mu = cf["var"], cf["mu"]
        P3 = cf["c3"] * (g(D, "D3", k) - TV["k3"][k]); P4 = cf["c4"] * (g(D, "g4row", k) - TV["k4"][k]); P = P3 + P4
        dt3 = tfit(g(D, "D3", k) - TV["k3"][k], 6 * mu * var)[0]; dt4 = tfit(g(D, "g4row", k) - TV["k4"][k], 12 * var * var)[0]
        pred = cf["c3"] * 6 * mu * var * dt3 + cf["c4"] * 12 * var * var * dt4
        cs = float(P @ pred) / max(np.linalg.norm(P) * np.linalg.norm(pred), 1e-300)
        line.append(f"{k}:{cs:+.2f}/{1 - float((P - pred) @ (P - pred)) / max(float(P @ P), 1e-300):+.2f}/{np.linalg.norm(pred) / max(np.linalg.norm(P), 1e-300):.2f}")
    print(f"    run {r} ({dumps[r]}): " + " ".join(line), flush=True)

print("  (4) null companions of the kappa4 error trace core dG (covariance metric at the true Sigma): readout-weighted cosines"
      "\n      of the actual errors with the companions [var: cv dvar vs cv dG_aa/12 | D3: c3 dk3 vs c3 mu dG_aa/4 |"
      "\n      C: Phi_a Phi_b dC vs Phi_a Phi_b dG_ab/12 | D21: rho_a Phi_b dD21/2 vs rho_a Phi_b (mu.dG)_aab/8]; and norm ratios")
G0 = {}
for r, D in enumerate(runs):
    out = []
    for k in range(1, L):
        if g(D, "g4row", k) is None or g(D, "wk431", k) is None:
            continue
        cf = coef(k); var, mu = cf["var"], cf["mu"]
        C = offd(sym(TM["cov"][k]))
        Gt, _, _ = fit_cov(TV["k4"][k].astype(np.float64), offd(TM["K22"][k]), offd(TM["K31"][k]), var, C)
        Gc, _, _ = fit_cov(g(D, "g4row", k), offd(g(D, "wk4m", k)), offd(g(D, "wk431", k).T), var, C)
        dG = Gc - Gt
        if r == 0:
            G0[k] = Gc
        dvar = g(D, "var", k) - var; dC = offd(sym(g(D, "C_off", k))) - C
        dk3 = g(D, "D3", k) - TV["k3"][k]; dD21 = offd(g(D, "D21", k)) - offd(TM["D21"][k])
        Pa = np.outer(cf["Ph"], cf["Ph"]); Ra = np.outer(cf["rho"], cf["Ph"])
        comp = {"var": (cf["cv"] * dvar, cf["cv"] * np.diag(dG) / 12),
                "D3": (cf["c3"] * dk3, cf["c3"] * mu * np.diag(dG) / 4),
                "C": (Pa * dC, Pa * offd(dG) / 12),
                "D21": (Ra * dD21 / 2, Ra * offd((2 * mu[:, None] * dG + mu[None, :] * np.diag(dG)[:, None]) / 3) / 8)}
        cell = []
        for nm, (x, y) in comp.items():
            cs = float(np.sum(x * y)) / max(math.sqrt(float(np.sum(x * x)) * float(np.sum(y * y))), 1e-300)
            cell.append(f"{nm} {cs:+.2f} ({math.sqrt(float(np.sum(y * y)) / max(float(np.sum(x * x)), 1e-300)):.2f})")
        if r > 0 and k in G0:
            dGr = Gc - G0[k]; Dk3 = g(D, "D3", k) - g(runs[0], "D3", k)
            x, y = cf["c3"] * Dk3, cf["c3"] * mu * np.diag(dGr) / 4
            cs = float(x @ y) / max(np.linalg.norm(x) * np.linalg.norm(y), 1e-300)
            cell.append(f"| change: D3 vs companion of dG change {cs:+.2f} ({np.linalg.norm(y) / max(np.linalg.norm(x), 1e-300):.2f})")
        out.append(f"      {k:2d}: " + "  ".join(cell))
    print(f"    run {r} ({dumps[r]}):\n" + "\n".join(out), flush=True)
