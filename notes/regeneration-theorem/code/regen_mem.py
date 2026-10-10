# Note XLIX section 8: the memory ladder for the (3,1) slice. The newborn leading structures of every layer m <= L,
# each evaluated with W_(m+1) replaced by the first-jet transported Wt_m = T_(L+1) ... T_(m+2) W_(m+1), summed.
#   python regen_mem.py NET L MCDIR WDIR
import sys, os, re, math, numpy as np
from scipy.special import ndtr
net, L, MC, WD = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4]
SLICE = os.environ.get("SLICE", "31")
NB = os.environ.get("NB", "1") == "1"    # genuinely newborn only: no kappa4 hyperedges, site factors without kappa4
NPAIRS = int(os.environ.get("NPAIRS", "8000"))
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "regen4_fast.py")).read()
g = {"np": np, "math": math, "re": re, "ndtr": ndtr, "SQ2PI": math.sqrt(2 * math.pi), "EMAX": 4,
     "EXACT_FALLBACK": False}
f64 = lambda a: np.asarray(a, dtype=np.float64)
T = {h: np.load(f"{MC}/mc2_off{net}_{h}.npz") for h in ("full", "h0", "h1")}
Wall = np.load(f"{WD}/W_off{net}.npy", mmap_mode="r")
n = Wall.shape[1]; g["n"] = n
exec(src[src.index("def he("):src.index("def r2(")], g)
full = {k: f64(T["full"][k]) for k in ("mu", "var", "k3", "k4")}
def layer_mats(m):
    C = f64(T["full"]["cov"][m]); Co = 0.5 * (C + C.T); np.fill_diagonal(Co, 0.0)
    d = {"Co": Co}
    for k in ("D21", "K31"):
        X = f64(T["full"][k][m]).copy(); np.fill_diagonal(X, 0.0); d[k] = X
    X = f64(T["full"]["K22"][m]); X = 0.5 * (X + X.T); np.fill_diagonal(X, 0.0); d["K22"] = X
    return d
Phi = {m: g["site_factors"](full["mu"][m], full["var"][m], full["k3"][m], full["k4"][m], True)[(1, 1)]
       for m in range(L + 1)}          # d kappa_1 = P(z > 0) under the true marginal
# transported matrices Wt_m, from m = L down to 0
Wt = {L: f64(Wall[L + 1])}
P = f64(Wall[L + 1])                    # running product T_(L+1) ... ; starts as W_(L+1)
for m in range(L - 1, -1, -1):
    P = (P * Phi[m + 1][None, :]) @ f64(Wall[m + 1])   # P <- P diag(Phi^(m+1)) W_(m+1)
    Wt[m] = P
mu1, v1 = full["mu"][L + 1], full["var"][L + 1]
s1 = np.sqrt(v1); a1 = mu1 / s1; act = a1 > -2.5
mask = (~np.eye(n, dtype=bool)) & act[:, None] & act[None, :]
ph1 = np.exp(-0.5 * a1 * a1) / math.sqrt(2 * math.pi)
rng = np.random.default_rng(1000 * net + L)
ia, ja = np.nonzero(mask); pick = rng.choice(len(ia), size=min(NPAIRS, len(ia)), replace=False); PI, PJ = ia[pick], ja[pick]
if SLICE == "31":
    rows, outdiag, pairs = "iiij", False, None
    sel = lambda X: X[mask]
    w = ((a1 * ph1 / s1**2)[:, None] * ndtr(a1)[None, :])[mask]
    tv = [f64(T[h]["K31"][L + 1])[mask] for h in ("h0", "h1", "full")]
    dil = (v1[:, None] * f64(T["full"]["cov"][L + 1]))[mask]
elif SLICE == "22":
    rows, outdiag, pairs = "iijj", False, (PI, PJ)
    sel = lambda X: X
    w = (ph1 / s1)[PI] * (ph1 / s1)[PJ]
    tv = []
    for h in ("h0", "h1", "full"):
        X = f64(T[h]["K22"][L + 1]); X = 0.5 * (X + X.T); tv.append(X[PI, PJ])
    dil = (v1[PI] * v1[PJ])
else:
    rows, outdiag, pairs = "iiii", True, None
    sel = lambda X: X[act]
    w = np.ones(int(act.sum()))
    tv = [f64(T[h]["k4"][L + 1])[act] for h in ("h0", "h1", "full")]
    dil = (v1 * v1)[act]
g["mats_W"] = None
def r2v(y):
    ip = lambda A, B: float(np.sum(A * B * w * w))
    return 1.0 - ip(tv[0] - y, tv[1] - y) / ip(tv[0], tv[1]), ip(y, tv[2]) / ip(y, y)
print(f"=== network {net}, memory ladder into the {SLICE} slice of layer {L + 1} (genuinely newborn only: {NB})", flush=True)
cum = 0.0
structs = [st for st in g["structures"](rows, True) if st[6] >= 1 and not (NB and st[4] is not None and st[4][0] in ("K22", "K31"))]
for m in range(L, -1, -1):
    F = g["site_factors"](full["mu"][m], full["var"][m], full["k3"][m], 0 * full["k4"][m] if NB else full["k4"][m], True)
    mats = layer_mats(m); mats["W"] = Wt[m]
    Ym = 0.0
    g["mats_W"] = Wt[m]
    for st in structs:
        Ym = Ym + g["evaluate"](rows, st, F, mats, outdiag, pairs=pairs)
    ym = sel(Ym)
    cum = cum + ym
    Rm, sm = r2v(ym); Rc, sc = r2v(cum)
    ipw = lambda A, B: float(np.sum(A * B * w * w))
    Rd, _ = r2v(cum + dil * (ipw(dil, tv[2] - cum) / ipw(dil, dil)))
    print(f"  source layer {m:2d}: alone R2 {100 * Rm:6.1f}% scale {sm:+.3f} | ladder from {L} down to {m}: "
          f"R2 {100 * Rc:6.1f}% scale {sc:+.3f} | + best-scaled dilation R2 {100 * Rd:6.1f}%", flush=True)

np.savez(os.environ.get("SAVE", f"ladder_{SLICE}_net{net}_L{L}.npz"), cum=cum, t0=tv[0], t1=tv[1], tf=tv[2], w=w,
         mask=(mask if SLICE == "31" else np.zeros(1)))
