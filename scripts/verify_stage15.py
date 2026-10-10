"""Stage 15 checks: deformed products (Ebrahimi-Fard--Patras--Tapia--Zambotti), row orthogonality, the error cocycle on an
official network, the Clifford/arcsine lift of the gate gas, and the quantum-calibration arithmetic.

    python scripts/verify_stage15.py [DATA_DIR]      (DATA_DIR holds W_off0.npy, truth_off0.npz; section 2 skipped if absent)

Writes nothing; the transcript is notes/stage15/verify_stage15.txt.
"""
import sys, os, itertools, numpy as np
from math import factorial, comb
from scipy.special import ndtr
from scipy import integrate

SQ2PI = np.sqrt(2 * np.pi)
phi = lambda x: np.exp(-0.5 * np.asarray(x, float) ** 2) / SQ2PI
rng = np.random.default_rng(15)


def hdr(s):
    print("\n" + "=" * 100 + "\n" + s + "\n" + "=" * 100, flush=True)


# ------------------------------------------------------------------------------------------------------------------
# 2-D formal power series helpers (coefficient arrays c[a, b] of s^a t^b / (a! b!), truncated at total degree D)
# ------------------------------------------------------------------------------------------------------------------
def egf_to_ord(c):
    a = np.arange(c.shape[0]); f = np.array([factorial(int(k)) for k in a], float)
    return c / np.outer(f, f)


def ord_to_egf(c):
    a = np.arange(c.shape[0]); f = np.array([factorial(int(k)) for k in a], float)
    return c * np.outer(f, f)


def series_mul(x, y, D):
    out = np.zeros_like(x)
    for a in range(D + 1):
        for b in range(D + 1 - a):
            out[a, b] = np.sum(x[:a + 1, :b + 1] * y[a::-1, b::-1][: a + 1, : b + 1])
    return out


def series_log(M, D):
    """log of an ordinary 2-D series with M[0,0] = 1, via homogeneous degree recursion n K_n = n M_n - sum k K_k M_{n-k}."""
    K = np.zeros_like(M)
    deg = np.add.outer(np.arange(M.shape[0]), np.arange(M.shape[1]))
    hom = lambda X, n: np.where(deg == n, X, 0.0)
    for n in range(1, D + 1):
        acc = n * hom(M, n)
        for k in range(1, n):
            acc -= k * hom(series_mul(hom(K, k), hom(M, n - k), D), n)
        K += acc / n
    return K


def series_exp(K, D):
    """exp of an ordinary 2-D series with K[0,0] = 0: n E_n = sum_k k K_k E_{n-k}."""
    E = np.zeros_like(K); E[0, 0] = 1.0
    deg = np.add.outer(np.arange(K.shape[0]), np.arange(K.shape[1]))
    hom = lambda X, n: np.where(deg == n, X, 0.0)
    for n in range(1, D + 1):
        acc = np.zeros_like(K)
        for k in range(1, n + 1):
            acc += k * hom(series_mul(hom(K, k), hom(E, n - k), D), n)
        E += acc / n
    return E


# ==================================================================================================================
hdr("1. Deformed products: Appell (Wick) polynomials of a non-Gaussian pair, the connected-cumulant theorem, and the\n"
    "   wall-jet expansion of Cov(relu, relu)")
# X = h1(g1), Y = h2(g2), (g1, g2) standard bivariate normal with correlation rho0; h_i(g) = g + a_i (g^2 - 1)
rho0, a1, a2, c, d = 0.6, 0.15, -0.10, 0.3, -0.2
NQ = 48
gx, gw = np.polynomial.hermite_e.hermegauss(NQ); gw = gw / gw.sum()
G1 = gx[:, None] + 0 * gx[None, :]; U = gx[None, :] + 0 * gx[:, None]
G2 = rho0 * G1 + np.sqrt(1 - rho0 ** 2) * U; Wq = np.outer(gw, gw)
h1 = lambda g: g + a1 * (g * g - 1); h2 = lambda g: g + a2 * (g * g - 1)
X, Y = h1(G1), h2(G2)
D = 16
mom = np.array([[np.sum(Wq * X ** p * Y ** q) if p + q <= D else 0.0 for q in range(D + 1)] for p in range(D + 1)])
Kc = ord_to_egf(series_log(egf_to_ord(mom), D))           # joint cumulants kappa_{ab}
print(f"joint cumulants: k10 {Kc[1,0]:+.2e} k20 {Kc[2,0]:.4f} k11 {Kc[1,1]:.4f} k30 {Kc[3,0]:+.4f} k21 {Kc[2,1]:+.4f}"
      f" k12 {Kc[1,2]:+.4f} k22 {Kc[2,2]:+.4f} k31 {Kc[3,1]:+.4f}")


def appell(mvec, P):
    """Appell polynomials A_0..A_P of a law with moments mvec: A = mu^{-1} * id (EFPTZ Thm 1.1), coefficient rows."""
    inv = np.zeros(P + 1); inv[0] = 1.0                    # mu^{-1} under binomial convolution
    for k in range(1, P + 1):
        inv[k] = -sum(comb(k, j) * mvec[j] * inv[k - j] for j in range(1, k + 1))
    A = np.zeros((P + 1, P + 1))
    for p in range(P + 1):
        for k in range(p + 1):
            A[p, k] = comb(p, k) * inv[p - k]
    return A


P = 7
AX, AY = appell(mom[:, 0], P), appell(mom[0, :], P)
evalp = lambda A, x: sum(A[k] * x ** k for k in range(A.shape[0]))
# connected theorem: E[A_p(X) A_q(Y)] = p! q! [s^p t^q] exp(sum_{a,b>=1} k_ab s^a t^b / a! b!)
Kmix = Kc.copy(); Kmix[0, :] = 0; Kmix[:, 0] = 0
Emix = ord_to_egf(series_exp(egf_to_ord(Kmix), D))
worst = 0.0
for p in range(1, P + 1):
    for q in range(1, P + 1):
        if p + q > D: continue
        lhs = np.sum(Wq * evalp(AX[p], X) * evalp(AY[q], Y)); worst = max(worst, abs(lhs - Emix[p, q]) / max(1, abs(lhs)))
print(f"E[A_p(X)A_q(Y)] vs connected-cumulant exponential, p,q <= {P}: max rel. difference {worst:.1e}")
print(f"  e.g. E[A_1A_1] {Emix[1,1]:.4f}=k11, E[A_2A_1] {Emix[2,1]:+.4f}=k21, E[A_2A_2] {Emix[2,2]:.4f}"
      f" (Mehler would give 2 k11^2 = {2*Kc[1,1]**2:.4f}; the rest is k22 + ...)")
# polynomial f, g: Cov(f(X), g(Y)) = sum_{p,q>=1} E f^(p) E g^(q) E[A_p A_q] / p! q!  (exact)
fX = X ** 3 - 2 * X; gY = Y ** 4 + Y
lhs = np.sum(Wq * fX * gY) - np.sum(Wq * fX) * np.sum(Wq * gY)
dfx = [np.sum(Wq * X ** 3 - 2 * X * Wq), np.sum(Wq * (3 * X ** 2 - 2)), np.sum(Wq * 6 * X), 6.0]
dgy = [np.sum(Wq * (Y ** 4 + Y)), np.sum(Wq * (4 * Y ** 3 + 1)), np.sum(Wq * 12 * Y ** 2), np.sum(Wq * 24 * Y), 24.0]
rhs = sum(dfx[p] * dgy[q] * Emix[p, q] / (factorial(p) * factorial(q)) for p in range(1, 4) for q in range(1, 5))
print(f"Cov(X^3-2X, Y^4+Y): quadrature {lhs:.10f}, Appell/connected expansion {rhs:.10f}")

# ReLU: wall jets J_1 = P(X>c), J_p = (-1)^p p_X^(p-2)(c), from the exact density of X = h1(g1)
import sympy as sp
xs = sp.symbols('x')
def density_jets(a, cc, P):
    g = (-1 + sp.sqrt(1 + 4 * a * (xs + a))) / (2 * a)       # inverse of g + a(g^2-1) on the main branch
    p = sp.exp(-g ** 2 / 2) / sp.sqrt(2 * sp.pi) / (1 + 2 * a * g)
    gstar = float(g.subs(xs, cc))
    J = [None, float(ndtr(-gstar))]
    for k in range(2, P + 1):
        J.append(float((-1) ** k * sp.diff(p, xs, k - 2).subs(xs, cc)))
    return J, gstar
JX, gsx = density_jets(a1, c, P); JY, gsy = density_jets(a2, d, P)
# exact Cov(relu(X-c), relu(Y-d)) by the analytic inner integral over u
beta = np.sqrt(1 - rho0 ** 2)
def inner(g1):
    ustar = (gsy - rho0 * g1) / beta; Q = ndtr(-ustar); f = phi(ustar)
    m0, m1, m2 = Q, f, ustar * f + Q                          # int_{u*} u^k phi(u) du, k = 0,1,2
    b0 = rho0 * g1; # h2(b0 + beta u) - d = (b0 + beta u) + a2((b0 + beta u)^2 - 1) - d
    c0 = b0 + a2 * (b0 * b0 - 1) - d; c1 = beta + 2 * a2 * b0 * beta; c2 = a2 * beta * beta
    return c0 * m0 + c1 * m1 + c2 * m2
EXY = integrate.quad(lambda g: (h1(g) - c) * inner(g) * phi(g), gsx, np.inf, epsabs=1e-14, epsrel=1e-13, limit=400)[0]
EX = integrate.quad(lambda g: (h1(g) - c) * phi(g), gsx, np.inf, epsabs=1e-14)[0]
EY = integrate.quad(lambda g: (h2(g) - d) * phi(g), gsy, np.inf, epsabs=1e-14)[0]
cov_true = EXY - EX * EY
# Gaussian closure: moment-matched bivariate normal, Mehler with Gaussian jets
mx, my, sx, sy = Kc[1, 0], Kc[0, 1], np.sqrt(Kc[2, 0]), np.sqrt(Kc[0, 2]); rr = Kc[1, 1] / (sx * sy)
ax_, ay_ = (mx - c) / sx, (my - d) / sy
def gjet(al, s, k):
    if k == 1: return ndtr(al)
    return phi(al) * np.polynomial.hermite_e.hermeval(-al, np.eye(k - 1)[k - 2]) / s ** (k - 1)   # s^{1-k} He_{k-2}(-a) phi(a)
cov_gauss = sum((sx * sy) ** p * gjet(ax_, sx, p) * gjet(ay_, sy, p) * rr ** p / factorial(p) for p in range(1, 40))
print(f"\nCov(relu(X-c), relu(Y-d)): exact {cov_true:.8f}; Gaussian closure {cov_gauss:.8f} (error {cov_gauss-cov_true:+.2e})")
print(" wall jets of X at c:", " ".join(f"{j:+.4f}" for j in JX[1:]))
for K in range(1, P + 1):
    tj_meh = sum(JX[p] * JY[p] * Kc[1, 1] ** p / factorial(p) for p in range(1, K + 1))
    tj_full = sum(JX[p] * JY[q] * Emix[p, q] / (factorial(p) * factorial(q)) for p in range(1, K + 1) for q in range(1, K + 1))
    print(f"  order <= {K}: true jets + Mehler (k11 only) {tj_meh:.8f} ({tj_meh-cov_true:+.1e});"
          f"  true jets + all connected cumulants {tj_full:.8f} ({tj_full-cov_true:+.1e})")


# ==================================================================================================================
hdr("2. Row orthogonality (Mehler for Gaussian rows) and the error cocycle on official network 0")
# 2a. exact row-Mehler: for fixed a, b and w ~ N(0, I/n): E[He_p(w.a/|a|/s) He_q(w.b/|b|/s)] = delta_pq p! rho^p
n0 = 64; a = np.abs(rng.standard_normal(n0)); b = np.abs(a + 0.8 * rng.standard_normal(n0))
ah, bh = a / np.linalg.norm(a), b / np.linalg.norm(b); rho = ah @ bh
Wr = rng.standard_normal((400000, n0)); u, v = Wr @ ah, Wr @ bh
He = lambda k, x: np.polynomial.hermite_e.hermeval(x, np.eye(k + 1)[k])
prods = [[He(p, u) * He(q, v) for q in range(4)] for p in range(4)]
tab = np.array([[x.mean() for x in row] for row in prods]); se = np.array([[x.std() / np.sqrt(len(x)) for x in row] for row in prods])
pred = np.diag([factorial(p) * rho ** p for p in range(4)])
print(f"rho = {rho:.4f}; E[He_p He_q] - delta_pq p! rho^p over p,q<=3: max |z| = {np.max((np.abs(tab - pred) / np.where(se > 0, se, np.inf))):.2f} (MC, 4e5 rows)")
# the radial obstruction: normalizing by a fixed scale instead of |a| leaks degree 2 into degree 0
cfix = 0.9 * np.linalg.norm(a)
print(f"radial leak: E[He_2(w.a/c)] = |a|^2/c^2 - 1 = {np.linalg.norm(a)**2/cfix**2 - 1:.4f}, MC {np.mean(He(2, Wr @ a / cfix)):.4f}")

D_DIR = sys.argv[1] if len(sys.argv) > 1 else "/tmp/claude-0/-home-user-thehonesttorus-github-io/9523ef6a-ff73-5ab2-8be3-275df2b5fcc6/scratchpad/data"
KP = os.path.join(os.path.dirname(D_DIR.rstrip("/")), "kprop3", "kprop3_0.npz")
if os.path.exists(os.path.join(D_DIR, "W_off0.npy")):
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
    from whest.relu_gauss import gauss_chain
    Wn = np.load(os.path.join(D_DIR, "W_off0.npy")).astype(np.float64); tr = np.load(os.path.join(D_DIR, "truth_off0.npz"))
    mt = tr["m"]; noise = tr["avg_variance"] / 1e9
    rec = {}; outG = gauss_chain(Wn, K=10, record=rec)
    ests = [("Gaussian closure", outG, np.array([rec[l]["alpha"] for l in range(16)]))]
    if os.path.exists(KP):
        Xk = np.load(KP); ests.append(("exact first-order chain (+ radial k4)", Xk["out"], Xk["m"] / np.sqrt(Xk["var"])))
    for name, out, R in ests:
        dm = out - mt; tot = np.mean(dm ** 2, 1); loc = np.zeros(16); tf = np.zeros(16); col = np.zeros(16); loc[0] = tot[0]
        for l in range(1, 16):
            trn = ndtr(R[l]) * (Wn[l] @ dm[l - 1]); loc[l] = np.mean((dm[l] - trn) ** 2); tf[l] = np.mean(trn ** 2) / tot[l - 1]
            uu = mt[l - 1] / np.linalg.norm(mt[l - 1]); col[l] = (dm[l - 1] @ uu) ** 2 / (len(uu) * tot[l - 1])
        e = dm[15]; t1 = ndtr(R[15]) * (Wn[15] @ dm[14])
        gain = np.array([np.prod(tf[l + 1:]) for l in range(16)])
        print(f"\n{name}: final MSE {tot[15]:.3e} (truth MC noise ~{noise:.0e})")
        print(f"  first-jet share of the final MSE (inherited mean error x Phi(r)): {1 - np.mean((e - t1)**2)/tot[15]:.3f}")
        print(f"  cocycle sum rule: sum_l local_l x gain_l = {np.sum(loc*gain):.3e} vs final {tot[15]:.3e}")
        print("  local defect x downstream gain, layers 1..16:", " ".join(f"{x:.1e}" for x in loc * gain))
        print("  transport factor per layer:", " ".join(f"{x:.2f}" for x in tf[1:]))
        print("  collective share of the inherited error:", " ".join(f"{x:.2f}" for x in col[1:]))
        print(f"  share of the final MSE born in layers 11-16: {np.sum((loc*gain)[10:])/np.sum(loc*gain):.2f}")
        # the same cocycle with the collective (mean-direction) component projected out at every layer
        Pp = lambda x, l: x - (x @ (mt[l] / np.linalg.norm(mt[l]))) * (mt[l] / np.linalg.norm(mt[l]))
        dp = np.array([Pp(dm[l], l) for l in range(16)]); totp = np.mean(dp ** 2, 1); locp = np.zeros(16); tfp = np.zeros(16); locp[0] = totp[0]
        for l in range(1, 16):
            trn = Pp(ndtr(R[l]) * (Wn[l] @ dp[l - 1]), l); locp[l] = np.mean((dp[l] - trn) ** 2); tfp[l] = np.mean(trn ** 2) / totp[l - 1]
        gp = np.array([np.prod(tfp[l + 1:]) for l in range(16)])
        print(f"  collective mode removed: final {totp[15]:.3e} ({totp[15]/tot[15]:.2f} of MSE), sum rule {np.sum(locp*gp):.3e}")
else:
    print("(official network files not found; section 2b skipped)")


# ==================================================================================================================
hdr("3. The gate gas as the commutative shadow of a quasi-free (Clifford) state: arcsine pairs and star cumulants")
from scipy.stats import multivariate_normal
def sign4(R):
    """E[e1 e2 e3 e4] for e = sign(z), z ~ N(0, R), from the 16 orthant probabilities (Genz)."""
    tot = 0.0
    for sgn in itertools.product([1, -1], repeat=4):
        Ds = np.diag(sgn); pr = multivariate_normal.cdf(np.zeros(4), mean=np.zeros(4), cov=Ds @ R @ Ds, abseps=1e-10, releps=1e-10, maxpts=2000000)
        tot += np.prod(sgn) * pr   # P(sgn_i z_i < 0 for all i) = P(sign z = -sgn); prod(-sgn) = prod(sgn)
    return tot
base = rng.uniform(-1, 1, (4, 4)); base = (base + base.T) / 2
for lam in (0.1, 0.2, 0.3):
    R = np.eye(4) + lam * (base - np.diag(np.diag(base)))
    E2 = (2 / np.pi) * np.arcsin(R)
    pair = E2[0, 1] * E2[2, 3] + E2[0, 2] * E2[1, 3] + E2[0, 3] * E2[1, 2]
    k4 = sign4(R) - pair
    star = -(4 / np.pi ** 2) * sum(np.prod([R[i, j] for j in range(4) if j != i]) for i in range(4))
    print(f"lambda {lam}: E[e1e2e3e4] {sign4(R):+.6f}; arcsine pairing {pair:+.6f}; sign cumulant k4 {k4:+.3e};"
          f" leading star -(4/pi^2) sum_i prod_j rho_ij {star:+.3e} (ratio {k4/star:.3f})")


# ==================================================================================================================
hdr("4. Calibration against quantum mean estimation (arithmetic; n = 1024, L = 16, B = 2048 n^3)")
n, L = 1024, 16; Bud = 2048 * n ** 3; fwd = 2 * L * n * n; sig = np.sqrt(0.0787)
for mse in (1.56e-8, 1e-9):
    eps = np.sqrt(mse); mc = sig ** 2 / mse * fwd / Bud; qone = sig / eps * fwd / Bud; qall = np.sqrt(n) * sig / eps * fwd / Bud
    print(f"target raw MSE {mse:.2e}: classical Monte Carlo {mc:8.2f} B; quantum mean estimation, one output {qone:.3f} B,"
          f" all n outputs (sqrt(n) multi-observable) {qall:.2f} B; deterministic closures reach 1.56e-8 at ~0.2 B")


# ==================================================================================================================
hdr("5. Magic (non-Gaussianity) injected by a fold: m(r) = Var relu(r+xi) - Phi(r)^2, its integral and the layer law")
mfun = lambda r: (1 + r * r) * ndtr(r) + r * phi(r) - (r * ndtr(r) + phi(r)) ** 2 - ndtr(r) ** 2
Minf = integrate.quad(mfun, -np.inf, np.inf, epsabs=1e-13)[0]
print(f"m(0) = {mfun(0.0):.5f} (= 1/2 - 1/2pi - 1/4); M_inf = int m(r) dr = {Minf:.6f}")
byk = [integrate.quad(lambda r, k=k: np.polynomial.hermite_e.hermeval(r, np.eye(k - 1)[k - 2]) ** 2 * phi(r) ** 2 / factorial(k), -np.inf, np.inf)[0] for k in range(2, 9)]
print("  int E_k(r) dr, k = 2..8:", " ".join(f"{x:.4f}" for x in byk), f"(sum to 8: {sum(byk):.4f})")
for t in (1, 4, 12):
    ex = integrate.quad(lambda r: mfun(r) * phi(r / np.sqrt(t)) / np.sqrt(t), -np.inf, np.inf)[0]
    print(f"  field law r ~ N(0,{t}): E m(r) = {ex:.5f} vs M_inf/sqrt(2 pi (1+t)) = {Minf/np.sqrt(2*np.pi*(1+t)):.5f}")


# ==================================================================================================================
hdr("6. Last-layer birth channels of the exact chain: residual-variance profile fitted by squared jets")
if os.path.exists(os.path.join(D_DIR, "W_off0.npy")) and os.path.exists(KP):
    from scipy.optimize import nnls
    R15 = Xk["m"][15] / np.sqrt(Xk["var"][15]); s15 = np.sqrt(Xk["var"][15])
    e = Xk["out"][15] - mt[15]; t1 = ndtr(R15) * (Wn[15] @ (Xk["out"][14] - mt[14])); res = e - t1
    Bf = np.stack([s15 * phi(R15) * np.polynomial.hermite_e.hermeval(R15, np.eye(8)[k]) for k in range(6)] + [s15 * ndtr(R15), s15 * R15 * ndtr(R15)], 1)
    cf, *_ = np.linalg.lstsq(Bf, res, rcond=None); q = res - Bf @ cf
    print(f"last-layer remainder (after first jet): {np.mean(res**2):.3e}; field-explained {1-np.mean(q**2)/np.mean(res**2):.3f}")
    bins = np.quantile(R15, np.linspace(0, 1, 21)); idx = np.clip(np.digitize(R15, bins) - 1, 0, 19)
    rb = np.array([R15[idx == b].mean() for b in range(20)]); vb = np.array([np.mean(q[idx == b] ** 2) for b in range(20)])
    Jb = lambda r: np.stack([phi(r) ** 2, (r * phi(r)) ** 2, ((r * r - 1) * phi(r)) ** 2], 1)
    co, rn = nnls(Jb(rb), vb)
    ch = (Jb(R15) * co).mean(0)
    print(f"  squared-jet fit (phi^2: covariance, (r phi)^2: kappa3, ((r^2-1) phi)^2: kappa4): coefficients {co}, rel. residual {rn/np.linalg.norm(vb):.2f}")
    print(f"  channel energies over the actual fields: covariance {ch[0]:.2e}, kappa3 {ch[1]:.2e}, kappa4 {ch[2]:.2e}"
          f" -> shares {ch[0]/ch.sum():.2f} / {ch[1]/ch.sum():.2f} / {ch[2]/ch.sum():.2f}")
