# Ward-identity accounting of the two-gain ledger's row sums.  The scale tangent of a Gaussian reference N(mu, S) under
# y -> G y (E G^2 = 1, Var G^2 = t) decomposes exactly into reference drift (dmu = -mu/8, dS = mu mu^T / 4), the
# third-cumulant direction T3 = sm kappa3 slices and the fourth-cumulant direction T4 = sm kappa4 slices, with no fifth
# or higher component; positive homogeneity makes the total transport of the output cumulants exact (1 + O(g)).  The
# ledger's rows r33 + r34 and r43 + r44 are the T3 + T4 parts on pair slices and pair output classes.  Here the drift
# part is added by central finite differences of the Gaussian one loop, and for the kappa4 row the mixture's dropped
# output classes (k4classes_gauss.py) as well; the sums are compared with 1.
import numpy as np, sys, time
sys.path.insert(0, "../num12")
from closure import relu_coeffs, relu2_coeffs
from pairvar import Layer, y_cumulants
net = int(sys.argv[1]) if len(sys.argv) > 1 else 0; eps = 0.05
D = np.load(f"mc_cum_off{net}.npz"); L, n = D["s1y"].shape
Wc = np.load(f"../official/W_off{net}.npy").astype(np.float64)
Z = np.load(f"gac3_off{net}.npz"); rows, gm3, gm4 = Z["rows"], Z["gm3"], Z["gm4"]
def fit(x, y): return (x @ y) / (x @ x)
def coh_k3(W, k3h, C3):
    C3o = C3.copy(); np.fill_diagonal(C3o, 0.0); return (W**3) @ k3h + 3 * (((W * W) @ C3o) * W).sum(1)
def coh_k4(W, k4h, C4):
    C4o = C4.copy(); np.fill_diagonal(C4o, 0.0); P = W * W; return (W**4) @ k4h + 3 * ((P @ C4o) * P).sum(1)
def oneloop(mu, S, W, x3, x4):
    sig = np.sqrt(np.diag(S)); A = relu_coeffs(mu, sig, 19); B = relu2_coeffs(mu, sig, 19); lay = Layer(mu, S, A, B)
    return fit(x3, coh_k3(W, np.diag(lay.C3G).copy(), lay.C3G)), fit(x4, coh_k4(W, lay.k4G, lay.C4G))
print(f"net {net}: Ward accounting per unit gain t.  kappa3 row: r33 + r34 (T3 + T4, pair) + drift3 ;  kappa4 row: r43 + r44 + drift4 + dropped output classes of the mixture")
print("  l+1 | r33+r34 | drift3 | sum3 | alpha^2-weighted O(g) term (1/2 - a^2/4) g3 || r43+r44 | drift4 | dropped mix | sum4 | O(g) term 1.5 g3 a^2 + g4 (1 - a^2/2)")
t0 = time.time()
for r in rows:
    l = int(r[0]); r33, r34, r43, r44 = r[11], r[12], r[13], r[14]
    Y = y_cumulants(D, l); mu, S = Y["mu"], Y["S"]; Yn = y_cumulants(D, l + 1); x3 = 1.5 * Yn["mu"] * np.diag(Yn["S"]); x4 = 3 * np.diag(Yn["S"])**2
    W = Wc[l + 1]; dmu = -mu / 8; dS = np.outer(mu, mu) / 4
    f3p, f4p = oneloop(mu + eps * dmu, S + eps * dS, W, x3, x4); f3m, f4m = oneloop(mu - eps * dmu, S - eps * dS, W, x3, x4)
    d3 = (f3p - f3m) / (2 * eps); d4 = (f4p - f4m) / (2 * eps)
    try: G = np.load(f"k4gauss_off{net}_l{l}.npz"); mix = float((G["m211"] + G["m1111"]) / gm4[l]) if gm4[l] > 1e-4 else float("nan")
    except Exception: mix = float("nan")
    a2 = np.diag(Yn["S"]); al2 = Yn["mu"]**2 / a2; w3 = fit(x3, x3 * (0.5 - al2 / 4)) * gm3[l + 1]; w4 = fit(x4, x4 * (1.5 * gm3[l + 1] * al2 + gm4[l + 1] * (1 - al2 / 2)))
    print(f"  {l+1:3d} | {r33+r34:.3f} | {d3:+.3f} | {r33+r34+d3:.3f} | {w3:+.3f} || {r43+r44:.3f} | {d4:+.3f} | {mix:+.3f} | {r43+r44+d4+mix:.3f} | {w4:+.3f}   [{time.time()-t0:.0f}s]", flush=True)
