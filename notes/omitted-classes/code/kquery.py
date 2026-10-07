# Checkpoint K's degree-one question in the chain's own terms (note XXXVI section 3i).
#   python kquery.py NET MCPREFIX DUMP0 [DUMP1 ...]      (free-running dumps with KEEP_EXTRA=pk1v,K2v,K11)
# The mean of layer k reads each pre-activation marginal only through its support query h(w_i) = E relu(z_i), a
# degree-one functional (checkpoint K, theorem 4.1). The chain represents it by the Edgeworth map of the physical
# cumulants (est_v29 TERM_SPECS[(1,)]):
#   m = Gauss(mu, s) + k3/6 E relu''' + k4/24 E relu'''' + k3^2/72 E relu^(6) + k3 k4/144 E relu^(7),
#   E relu^(j)(z) = He_(j-2)(-a) phi(a) / s^(j-1)  (j >= 2) for z ~ N(mu, s^2), a = mu/s.
# (0) check: the map reproduces the chain's own pk1v from its dumped mu, var, D3, g4row.
# (a) truncation: H_k = E relu(z_k) - map(true mu, var, k3, k4) is the part of the true query the kappa<=4 state
#     cannot represent even when exact. Its output image under the full mean-covariance propagation (fentry.py's
#     reference maps) bounds what exact cumulants would leave. Noise-free norms from the Monte Carlo halves.
#     If the hidden non-Gaussianity carried skewed scale (checkpoint K section 5), H would cancel part of the kappa4 term:
#     the cross-half regression coefficient of H on T4 = c4 k4 would be negative.
# (b) source level: per unit, the chain's query error decomposes as N_m = P3 + P4 + R with P3 = c3 dk3, P4 = c4 dk4
#     (dk = chain - truth, first-order coefficients at the truth) and R the rest (which contains -H). Do P3 and P4
#     compensate inside one query, and how do y2 / KD move that?
import sys, math, numpy as np
net, pre, dumps = int(sys.argv[1]), sys.argv[2], sys.argv[3:]
Wcol = np.load(f"../official/W_off{net}.npy").astype(np.float64); L, n, _ = Wcol.shape
mt = np.load(f"../official/truth_off{net}.npz")["m"].astype(np.float64)
eye = np.eye(n, dtype=bool)
erf = np.vectorize(math.erf)


def offd(X):
    X = np.array(X, dtype=np.float64); X[eye] = 0.0; return X


def gauss(mu, var):
    s = np.sqrt(var); a = mu / s
    ph = np.exp(-a * a / 2) / math.sqrt(2 * math.pi); Ph = 0.5 * (1 + erf(a / math.sqrt(2)))
    return s, a, ph, Ph


def emap(mu, var, k3, k4, parts=False):
    s, a, ph, Ph = gauss(mu, var)
    G = mu * Ph + s * ph
    d3 = -a * ph / s ** 2                                    # E relu'''  = He1(-a) phi / s^2
    d4 = (a * a - 1) * ph / s ** 3                           # E relu'''' = He2(-a) phi / s^3
    d6 = (a ** 4 - 6 * a * a + 3) * ph / s ** 5              # He4(-a)
    d7 = -(a ** 5 - 10 * a ** 3 + 15 * a) * ph / s ** 6      # He5(-a)
    T3, T4, T33, T34 = k3 / 6 * d3, k4 / 24 * d4, k3 * k3 / 72 * d6, k3 * k4 / 144 * d7
    m = G + T3 + T4 + T33 + T34
    return (m, G, T3, T4, T33, T34) if parts else m


T = {t: np.load(f"{pre}_{t}.npz") for t in ("full", "h0", "h1")}
tr = {t: {k: T[t][k].astype(np.float64) for k in ("mu", "var", "k3", "k4", "mu_y")} for t in T}
F = tr["full"]
sig, al, ph, Ph = gauss(F["mu"], F["var"])
cv = ph / (2 * sig); rho = ph / sig; dPdv = -al * ph / (2 * F["var"])
c3 = -al * ph / (6 * F["var"]); c4 = (al * al - 1) * ph / (24 * F["var"] * sig)
COV = np.load(f"{pre}_full.npz")["cov"]
Cref = [offd(0.5 * (COV[k] + COV[k].T)) for k in range(L)]
Kref = [np.outer(Ph[k], Ph[k]) + np.outer(rho[k], rho[k]) * Cref[k] for k in range(L)]


def forward_m(src):
    """Output image of m-level sources src[k] (n-vectors) under the reference maps (as fentry.py)."""
    dm, dC = np.zeros(n), np.zeros((n, n))
    for k in range(L):
        W = Wcol[k]
        dmu = W @ dm if k > 0 else np.zeros(n)
        dCz = W @ dC @ W.T if k > 0 else np.zeros((n, n))
        dvar = np.diag(dCz).copy(); m = mt[k]
        bm = Ph[k] * dmu + cv[k] * dvar + src[k]
        if k < L - 1:
            dv = 2 * m * (1 - Ph[k]) * dmu + (Ph[k] - 2 * m * cv[k]) * dvar
            u = rho[k] * dmu + dPdv[k] * dvar
            bC = Kref[k] * offd(dCz) + offd(np.outer(u, Ph[k]) * Cref[k] + np.outer(Ph[k], u) * Cref[k])
            dC = bC + np.diag(dv)
        dm = bm
    return dm


def adjoint_m(w):
    """m-level adjoint a_k of <w, m_(L-1)> under the reference maps (with the covariance channel)."""
    a, Lam = w.copy(), np.zeros((n, n)); res = [None] * L
    for k in range(L - 1, -1, -1):
        res[k] = a
        lam = np.diag(Lam).copy(); Lo = offd(Lam); h = (Lo * Cref[k]) @ Ph[k]; m = mt[k]
        alpha = a * Ph[k] + lam * 2 * m * (1 - Ph[k]) + 2 * rho[k] * h
        beta = a * cv[k] + lam * (Ph[k] - 2 * m * cv[k]) + 2 * dPdv[k] * h
        if k > 0:
            a = Wcol[k].T @ alpha; Lam = Wcol[k].T @ (Lo * Kref[k] + np.diag(beta)) @ Wcol[k]
    return res


runs = [np.load(p) for p in dumps]
gk = lambda D, key, k: D[f"{key}_{k}"].astype(np.float64) if f"{key}_{k}" in D.files else None
es = [D[f"pk1v_{L - 1}"].astype(np.float64) - mt[-1] for D in runs]
E0 = float(es[0] @ es[0])
print(f"net {net}: support queries E relu(z_i) and the kappa<=4 Edgeworth map; baseline n MSE {E0:.4e}", flush=True)

# (0) the map reproduces the chain's own means
D = runs[0]; worst = 0.0
for k in range(1, L):
    mu_c = Wcol[k] @ gk(D, "pk1v", k - 1)
    pred = emap(mu_c, gk(D, "var", k), gk(D, "D3", k), gk(D, "g4row", k))
    worst = max(worst, np.linalg.norm(pred - gk(D, "pk1v", k)) / np.linalg.norm(gk(D, "pk1v", k) - mt[k]))
print(f"  (0) map vs the chain's pk1v from its own dumped state: max |diff| / |chain error| = {worst:.2e}", flush=True)

# (a) truncation remainder H of the true queries
H = {t: [tr[t]["mu_y"][k] - emap(tr[t]["mu"][k], tr[t]["var"][k], tr[t]["k3"][k], tr[t]["k4"][k]) for k in range(L)] for t in T}
PT = {t: [emap(tr[t]["mu"][k], tr[t]["var"][k], tr[t]["k3"][k], tr[t]["k4"][k], parts=True) for k in range(L)] for t in T}
nf = lambda x0, x1: float(x0 @ x1)                          # noise-free inner product from the two halves
print("  (a) per layer, rms over units (noise-free): Gauss-to-truth gap, the map's kappa terms, the remainder H; and the\n"
      "      cross-half coefficient of H on T4 = c4 k4 and on T3 (b < 0: H cancels part of the kappa4 term)")
print("    layer: |true - Gauss|  |T3|  |T4|  |T33+T34|  |H|  | H/T4 coef  H/T3 coef | cos(H,T4) cos(H,T3)")
r = lambda v: math.sqrt(max(v, 0.0) / n)
for k in range(L):
    p0, p1 = PT["h0"][k], PT["h1"][k]
    gap0, gap1 = tr["h0"]["mu_y"][k] - p0[1], tr["h1"]["mu_y"][k] - p1[1]
    H0, H1 = H["h0"][k], H["h1"][k]
    T30, T31, T40, T41 = p0[2], p1[2], p0[3], p1[3]
    S0, S1 = p0[4] + p0[5], p1[4] + p1[5]
    hh, t44, t33 = nf(H0, H1), nf(T40, T41), nf(T30, T31)
    h4 = 0.5 * (nf(H0, T41) + nf(H1, T40)); h3 = 0.5 * (nf(H0, T31) + nf(H1, T30))
    co = lambda a_, b_, c_: a_ / math.sqrt(b_ * c_) if b_ > 0 and c_ > 0 else float("nan")
    print(f"    {k:2d}: {r(nf(gap0, gap1)):.3e} {r(t33):.3e} {r(t44):.3e} {r(nf(S0, S1)):.3e} {r(hh):.3e} | "
          f"{h4 / t44 if t44 > 0 else float('nan'):+6.3f} {h3 / t33 if t33 > 0 else float('nan'):+6.3f} | "
          f"{co(h4, hh, t44):+.3f} {co(h3, hh, t33):+.3f}", flush=True)
img = {t: forward_m([-H[t][k] for k in range(L)]) for t in T}       # the error exact cumulants would leave
print(f"  (a) output image of -H (what exact kappa3, kappa4 would leave), % of the baseline n MSE:")
print(f"      |image|^2 noise-free {100 * nf(img['h0'], img['h1']) / E0:+.2f}  (full truth {100 * float(img['full'] @ img['full']) / E0:.2f}, "
      f"half-split noise {100 * abs(float(img['h0'] @ img['h0']) - float(img['h1'] @ img['h1'])) / 2 / E0:.2f})")
for r_, e in enumerate(es):
    sh = [100 * float(img[t] @ e) / float(e @ e) for t in ("full", "h0", "h1")]
    print(f"      run {r_} ({dumps[r_]}): share <image, e>/|e|^2 {sh[0]:+.2f} +- {abs(sh[1] - sh[2]) / 2:.2f}")
adj = adjoint_m(es[0])
print("      per layer share of the baseline error, <a_k, -H_k>/|e0|^2 (%): " +
      " ".join(f"{100 * float(adj[k] @ (-H['full'][k])) / E0:+.2f}" for k in range(L)), flush=True)

# (b) source-level decomposition of each run's query error
print("  (b) per run and layer: rms of P3 = c3 dk3, P4 = c4 dk4, their sum, R = N_m - P3 - P4 (noise-free); cos(P3, P4)")
for r_, D in enumerate(runs):
    print(f"    run {r_} ({dumps[r_]}):  layer: |P3| |P4| |P3+P4| |R| |N_m| | cos(P3,P4) cos(P3+P4,R)")
    acc = np.zeros(4)
    for k in range(1, L):
        Pm = {}
        for t in ("h0", "h1"):
            dmu = Wcol[k] @ (gk(D, "pk1v", k - 1) - mt[k - 1])
            dvar = gk(D, "var", k) - tr[t]["var"][k]
            Nm = gk(D, "pk1v", k) - mt[k] - (Ph[k] * dmu + cv[k] * dvar)
            P3 = c3[k] * (gk(D, "D3", k) - tr[t]["k3"][k]); P4 = c4[k] * (gk(D, "g4row", k) - tr[t]["k4"][k])
            Pm[t] = (P3, P4, P3 + P4, Nm - P3 - P4, Nm)
        v = [nf(Pm["h0"][j], Pm["h1"][j]) for j in range(5)]
        c34 = 0.5 * (nf(Pm["h0"][0], Pm["h1"][1]) + nf(Pm["h1"][0], Pm["h0"][1]))
        cSR = 0.5 * (nf(Pm["h0"][2], Pm["h1"][3]) + nf(Pm["h1"][2], Pm["h0"][3]))
        co = lambda a_, b_, c_: a_ / math.sqrt(b_ * c_) if b_ > 0 and c_ > 0 else float("nan")
        acc += np.array([v[0], v[1], c34, v[2]])
        print(f"      {k:2d}: {r(v[0]):.3e} {r(v[1]):.3e} {r(v[2]):.3e} {r(v[3]):.3e} {r(v[4]):.3e} | "
              f"{co(c34, v[0], v[1]):+.3f} {co(cSR, v[2], v[3]):+.3f}", flush=True)
    print(f"      all layers: |P3|^2 {acc[0]:.3e} |P4|^2 {acc[1]:.3e} 2<P3,P4> {2 * acc[2]:+.3e} |P3+P4|^2 {acc[3]:.3e}")
    if r_ > 0:
        # source-level analogue of the pair terms: does the change of P3 line up with the baseline's P4?
        line = []
        for k in range(1, L):
            dP3 = c3[k] * (gk(D, "D3", k) - gk(runs[0], "D3", k)); dP4 = c4[k] * (gk(D, "g4row", k) - gk(runs[0], "g4row", k))
            P4b = 0.5 * (c4[k] * (gk(runs[0], "g4row", k) - tr["h0"]["k4"][k]) + c4[k] * (gk(runs[0], "g4row", k) - tr["h1"]["k4"][k]))
            P3b = 0.5 * (c3[k] * (gk(runs[0], "D3", k) - tr["h0"]["k3"][k]) + c3[k] * (gk(runs[0], "D3", k) - tr["h1"]["k3"][k]))
            line.append(f"{k}:{float(dP3 @ P4b) / max(np.linalg.norm(dP3) * np.linalg.norm(P4b), 1e-300):+.2f}/"
                        f"{float(dP4 @ P3b) / max(np.linalg.norm(dP4) * np.linalg.norm(P3b), 1e-300):+.2f}")
        print("      cos(change of P3, baseline P4) / cos(change of P4, baseline P3) per layer: " + " ".join(line), flush=True)
