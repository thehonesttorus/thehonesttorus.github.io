# Validation of the regen4 engine on lower orders with known answers: covariance (legs ij), D21 (legs iij),
# kappa3 diagonal (legs iii) of z' = W relu(z), against the layer-(L+1) truth.
import sys
net, L, MC, WD = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4]
EMAXV = int(sys.argv[5]) if len(sys.argv) > 5 else 3
src = open(__file__.replace("validate_low.py", "regen4.py")).read()
head = src[:src.index("lay = {h:")]
sys.argv = ["regen4.py", str(net), str(L), MC, WD]
exec(head)
EMAX = EMAXV
lay = {k: f64(T["full"][k][L]) for k in ("mu", "var", "k3", "k4", "cov", "D21", "K22", "K31")}
nxt = {h: {k: f64(T[h][k][L + 1]) for k in ("mu", "var", "k3", "D21", "cov")} for h in ("full", "h0", "h1")}
Co = 0.5 * (lay["cov"] + lay["cov"].T); np.fill_diagonal(Co, 0.0)
mats = {"W": W, "Co": Co}
for k in ("D21", "K31"):
    X = lay[k].copy(); np.fill_diagonal(X, 0.0); mats[k] = X
X = 0.5 * (lay["K22"] + lay["K22"].T); np.fill_diagonal(X, 0.0); mats["K22"] = X
Ft = site_factors(lay["mu"], lay["var"], lay["k3"], lay["k4"], True)
a1 = nxt["full"]["mu"] / np.sqrt(nxt["full"]["var"]); act = a1 > -2.5
mask = off & act[:, None] & act[None, :]
for name, rows, tk, diag in (("cov offdiag", "ij", "cov", False), ("var diag", "ii", "var", True),
                             ("D21", "iij", "D21", False), ("k3 diag", "iii", "k3", True)):
    Tk = []
    for h in ("h0", "h1", "full"):
        X = nxt[h][tk].copy()
        if not diag:
            np.fill_diagonal(X, 0.0)
        Tk.append(X)
    msk = act if diag else mask
    tot = 0.0; by = {}
    for st in structures(rows, hyper=True):
        Y = evaluate(rows, st, Ft, mats, diag)
        if not diag:
            np.fill_diagonal(Y, 0.0)
        e = sum(st[3].values()); key = ("h:" + st[4][0] if st[4] else "C") + f" e={e}"
        by[key] = by.get(key, 0.0) + Y; tot = tot + Y
    ip = lambda A, B: float(np.sum((A * B)[msk]))
    def R(Y):
        return 1 - ip(Tk[0] - Y, Tk[1] - Y) / ip(Tk[0], Tk[1]), ip(Y, Tk[2]) / ip(Y, Y)
    r, s = R(tot)
    print(f"{name}: all R2 {100*r:.3f}% scale {s:+.4f}", flush=True)
    for k in sorted(by):
        r, s = R(by[k]); print(f"    {k:12s} alone R2 {100*r:8.3f}% scale {s:+.4f}", flush=True)
