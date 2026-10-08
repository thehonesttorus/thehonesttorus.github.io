# The readout at true inputs, all three channels (note XXXVI section 3j).
#   python cread.py NET MCPREFIX ORACLE_A ORACLE_B BASEDUMP
# ORACLE_A / ORACLE_B: all-oracle runs from the two independent halves of a second Monte Carlo set; each is an all-oracle run (oracle_one.py NET dump:D3+D21+G4+WK4M+K31+VAR+COFF <independent Monte Carlo> ...,
# KEEP_EXTRA=K11,pk1v,K2v): at every layer the chain reads the pre-activation slices, variance and covariance from an
# independent Monte Carlo set, so its post-activation mean, variance and covariance are its readout (Edgeworth mean map +
# gated pair program) of (nearly) true inputs. BASEDUMP: the production free-running dump (for the output error e).
# Readout defect at layer k, removing the response to the small input deviations with the reference linearisation B_k of
# fentry.py:  N^read_k = Delta X_k - B_k Delta Z_k,  Delta = oracle run - truth.
# Its parts (mean, var_y, C^y offdiag) are propagated by the full mean-covariance maps; their output images say how much of
# the current error the readout alone would leave (energy, noise-free from the two halves of the truth) and how much of the
# baseline error it carries (share <image, e>/|e|^2). Noise-free energies pair (ORACLE_A, truth h0) with (ORACLE_B, truth
# h1): all four Monte Carlo sets are independent. Section 3i measured the mean part analytically (H); this adds the
# second-moment readout (the pair program at true inputs).
import sys, math, numpy as np
net, pre, opa, opb, bp = int(sys.argv[1]), sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]
Wcol = np.load(f"../official/W_off{net}.npy").astype(np.float64); L, n, _ = Wcol.shape
mt = np.load(f"../official/truth_off{net}.npz")["m"].astype(np.float64)
eye = np.eye(n, dtype=bool); erf = np.vectorize(math.erf)


def offd(X):
    X = np.array(X, dtype=np.float64); X[eye] = 0.0; return X


def sym(X):
    X = np.asarray(X, dtype=np.float64); return 0.5 * (X + X.T)


T = {t: np.load(f"{pre}_{t}.npz") for t in ("full", "h0", "h1")}
F = T["full"]
mu, var = F["mu"].astype(np.float64), F["var"].astype(np.float64)
sig = np.sqrt(var); al = mu / sig; ph = np.exp(-al * al / 2) / math.sqrt(2 * math.pi); Ph = 0.5 * (1 + erf(al / math.sqrt(2)))
cv = ph / (2 * sig); rho = ph / sig; dPdv = -al * ph / (2 * var)
COV = F["cov"]
Cref = [offd(sym(COV[k])) for k in range(L)]
Kref = [np.outer(Ph[k], Ph[k]) + np.outer(rho[k], rho[k]) * Cref[k] for k in range(L)]
TR = {t: dict(var=T[t]["var"].astype(np.float64), vy=T[t]["var_y"].astype(np.float64), cov=T[t]["cov"], cy=T[t]["cov_y"]) for t in T}


def B(k, dmu, dCz):
    dvar = np.diag(dCz).copy(); m = mt[k]
    dm = Ph[k] * dmu + cv[k] * dvar
    if k == L - 1:
        return dm, None, None
    dv = 2 * m * (1 - Ph[k]) * dmu + (Ph[k] - 2 * m * cv[k]) * dvar
    u = rho[k] * dmu + dPdv[k] * dvar
    return dm, dv, Kref[k] * offd(dCz) + offd(np.outer(u, Ph[k]) * Cref[k] + np.outer(Ph[k], u) * Cref[k])


def forward(src):
    """src[k] = (sm, sv, sC) y-level sources (None allowed); returns the output image."""
    dm, dC = np.zeros(n), np.zeros((n, n))
    for k in range(L):
        W = Wcol[k]
        dmu = W @ dm if k > 0 else np.zeros(n)
        dCz = W @ dC @ W.T if k > 0 else np.zeros((n, n))
        bm, bv, bC = B(k, dmu, dCz)
        sm, sv, sC = src[k]
        if sm is not None: bm = bm + sm
        if k < L - 1:
            if sv is not None: bv = bv + sv
            if sC is not None: bC = bC + sC
            dC = bC + np.diag(bv)
        dm = bm
    return dm


OA, OB, Bd = np.load(opa), np.load(opb), np.load(bp)
PAIRS = {"full": (OA, "full"), "h0": (OA, "h0"), "h1": (OB, "h1")}
g = lambda D, key, k: D[f"{key}_{k}"].astype(np.float64) if f"{key}_{k}" in D.files else None
e = g(Bd, "pk1v", L - 1) - mt[-1]; E = float(e @ e)
print(f"net {net}: readout of (nearly) true inputs, all channels; oracle runs {opa}, {opb}; baseline n MSE {E:.4e}")
for O_ in (OA, OB):
    print(f"  oracle run output n MSE {float((g(O_, 'pk1v', L - 1) - mt[-1]) @ (g(O_, 'pk1v', L - 1) - mt[-1])):.4e}")
src = {t: [] for t in T}; diag = []
for k in range(L):
    W = Wcol[k]
    for t in T:
        tt = TR[t]; O = PAIRS[t][0]
        Czt = offd(sym(tt["cov"][k])) + np.diag(tt["var"][k])
        if k == 0:
            dmu = np.zeros(n)
        else:
            dmu = W @ (g(O, "pk1v", k - 1) - mt[k - 1])
        Cz = offd(sym(g(O, "C_off", k))) + np.diag(g(O, "var", k))
        dCz = Cz - Czt
        dm = g(O, "pk1v", k) - mt[k]
        bm, bv, bC = B(k, dmu, dCz)
        if k < L - 1:
            Cy = offd(sym(g(O, "K11", k))) + np.diag(g(O, "K2v", k))
            dCy = Cy - (offd(sym(tt["cy"][k])) + np.diag(tt["vy"][k]))
            sv, sC = np.diag(dCy) - bv, offd(dCy) - bC
            if t == "full":
                diag.append((k, np.linalg.norm(offd(dCy)) / np.linalg.norm(offd(sym(tt["cy"][k]))), np.linalg.norm(sC) / max(np.linalg.norm(offd(dCy)), 1e-300),
                             np.linalg.norm(bC) / max(np.linalg.norm(offd(dCy)), 1e-300)))
        else:
            sv = sC = None
        src[t].append((dm - bm, sv, sC))
print("  per layer: |dC^y_off|/|C^y_off| of the oracle run, |readout defect|/|dC^y|, |input-response part|/|dC^y|")
for k, a, b, c in diag:
    print(f"    {k:2d}: {a:.4f} {b:.3f} {c:.3f}")
parts = {"mean": (0,), "var_y": (1,), "C_off": (2,), "all": (0, 1, 2)}
img = {}
for name, idx in parts.items():
    for t in T:
        s = [tuple(x[j] if j in idx else None for j in range(3)) for x in src[t]]
        img[(name, t)] = forward(s)
print("  output images of the readout defect, % of the baseline n MSE: energy noise-free (half-split noise) | share <image, e>/|e|^2")
for name in parts:
    i0, i1, f_ = img[(name, "h0")], img[(name, "h1")], img[(name, "full")]
    en = float(i0 @ i1); nz = abs(float(i0 @ i0) - float(i1 @ i1)) / 2
    sh = [float(img[(name, t)] @ e) for t in ("full", "h0", "h1")]
    print(f"    {name:6s}: {100 * en / E:+7.2f} ({100 * nz / E:.2f}) | {100 * sh[0] / E:+6.2f} +- {100 * abs(sh[1] - sh[2]) / 2 / E:.2f}", flush=True)
