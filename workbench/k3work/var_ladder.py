# Is the chain's per-layer variance injection perturbative? (note XLIII section 3)
#   python var_ladder.py LOCDIR MC2DIR NET "LAYERS"      (LAYERS = target layers t, e.g. "6,9,12,13,14,15")
# For source s = t - 1 and the true law of h_s (16M-input Monte Carlo, full and halves A, B):
#   g   = var_true_t - diag W_t KG(true mu_s, C_s) W_t^T        (true defect of the Gaussian closure)
#   E3  = diag W_t D3 W_t^T,  E4 = diag W_t D4 W_t^T             (first-order bivariate Edgeworth corrections of the
#         post-activation covariance from the true cumulants: order 3 = kappa(a,a,b), kappa3(a); order 4 =
#         kappa(a,a,b,b), kappa(a,a,a,b), kappa4(a); exact bivariate-normal coefficients, Hermite series in rho)
#   x   = var_chain_t - diag W_t KG(chain state_s) W_t^T         (the chain's own correction)
# Uses E_G[relu^(p)(h_a) relu^(q)(h_b)] = s_a^(1-p) s_b^(1-q) sum_m c_(m+p)(al_a) c_(m+q)(al_b) rho^m / m!, with
# c_m(al) = E[relu(al + z) He_m(z)]. Saves per-neuron arrays to var_ladder_{net}.npz for structure regressions.
import sys, numpy as np
from math import factorial
from scipy.special import ndtr

loc, mcd, net = sys.argv[1], sys.argv[2], int(sys.argv[3])
layers = [int(v) for v in sys.argv[4].split(",")]
M = 40
SQ2PI = np.sqrt(2.0 * np.pi)
FACT = np.array([float(factorial(k)) for k in range(M + 6)])


def ccoef(al, mmax):
    """c[m] = E[relu(al + z) He_m(z)] for m = 0..mmax (unnormalised)."""
    phi = np.exp(-0.5 * al * al) / SQ2PI; Phi = ndtr(al)
    c = np.zeros((mmax + 1, al.size))
    c[0] = al * Phi + phi; c[1] = Phi
    He_prev, He = np.zeros_like(al), np.ones_like(al)      # He_{-1}, He_0
    for m in range(2, mmax + 1):
        k = m - 2
        if k >= 1:
            He_prev, He = He, al * He - (k - 1.0) * He_prev
        c[m] = (-1.0) ** m * phi * He
    return c


def series(ca, cb, R, p, q, m0):
    """S_ab = sum_{m >= m0} ca[m+p]_a cb[m+q]_b R_ab^m / m!   (n x n)."""
    S = np.zeros_like(R); P = np.ones_like(R) if m0 == 0 else R.copy()
    for m in range(m0, M + 1):
        if m > m0:
            P = P * R
        S += np.outer(ca[m + p], cb[m + q]) * P / FACT[m]
    return S


def relu_var(mu, v):
    s = np.sqrt(v); al = mu / s; Phi = ndtr(al); phi = np.exp(-0.5 * al * al) / SQ2PI
    m1 = s * (al * Phi + phi); m2 = v * ((al * al + 1.0) * Phi + al * phi)
    return m2 - m1 * m1, m1


def KG(mu, C):
    v = np.diag(C).copy(); s = np.sqrt(v); al = mu / s
    R = C / np.outer(s, s); np.fill_diagonal(R, 0.0)
    c = ccoef(al, M + 1)
    K = np.outer(s, s) * series(c, c, R, 0, 0, 1)
    K[np.diag_indices_from(K)] = relu_var(mu, v)[0]
    return K


def edgeworth(mu, C, k3, k4, D21, K22, K31):
    """first-order Edgeworth corrections (D3, D4) to Cov(relu(h)) from the given cumulants."""
    v = np.diag(C).copy(); s = np.sqrt(v); al = mu / s
    R = C / np.outer(s, s); np.fill_diagonal(R, 0.0)
    c = ccoef(al, M + 5)
    sa = s[:, None]; sb = s[None, :]
    # order 3
    D3 = 0.5 * D21 / sa * series(c, c, R, 2, 1, 0)                     # kappa(a,a,b): (p,q) = (2,1)
    D3 = D3 + D3.T                                                      # (1,2) by symmetry of the roles
    T30 = (k3[:, None] / 6.0) * sa ** -2 * sb * series(c, c, R, 3, 0, 1)
    D3 += T30 + T30.T
    # order 4
    D4 = 0.25 * K22 / (sa * sb) * series(c, c, R, 2, 2, 0)              # kappa(a,a,b,b)
    T31 = (K31 / 6.0) * sa ** -2 * series(c, c, R, 3, 1, 0)             # kappa(a,a,a,b)
    D4 += T31 + T31.T
    T40 = (k4[:, None] / 24.0) * sa ** -3 * sb * series(c, c, R, 4, 0, 1)
    D4 += T40 + T40.T
    # diagonal: Var(relu(h_a)) corrections
    _, m1 = relu_var(mu, v)
    d3 = k3 / 6.0 * 2.0 * c[2] / s - 2.0 * m1 * k3 / 6.0 * c[3] / s ** 2
    d4 = k4 / 24.0 * 2.0 * c[3] / s ** 2 - 2.0 * m1 * k4 / 24.0 * c[4] / s ** 3
    D3[np.diag_indices_from(D3)] = d3; D4[np.diag_indices_from(D4)] = d4
    return D3, D4


def qdiag(W, K):
    return np.einsum("ij,ij->i", W @ K, W)


cd = np.load(f"{loc}/chaindump_{net}.npz")
W = np.load(f"{loc}/W_off{net}.npy").astype(np.float64)
T = {h: np.load(f"{mcd}/mc2_off{net}_{h}.npz") for h in ("full", "h0", "h1")}
out = cd["out"]
save = {}
print(f"=== network {net}: per-neuron quantities / true variance at the target layer; 'nf' = noise-free rms (halves)", flush=True)
for t in layers:
    s_ = t - 1
    mu_c = W[s_] @ out[s_ - 1]
    Cc = cd["C"][s_].astype(np.float64); Cc = 0.5 * (Cc + Cc.T); np.fill_diagonal(Cc, cd["var"][s_])
    x = cd["var"][t] - qdiag(W[t], KG(mu_c, Cc))
    res = {}
    for h in ("full", "h0", "h1"):
        F = T[h]
        mu = F["mu"][s_].astype(np.float64); C = F["cov"][s_].astype(np.float64); C = 0.5 * (C + C.T)
        g = F["var"][t].astype(np.float64) - qdiag(W[t], KG(mu, C))
        D21 = F["D21"][s_].astype(np.float64); K22 = F["K22"][s_].astype(np.float64)
        K22 = 0.5 * (K22 + K22.T); K31 = F["K31"][s_].astype(np.float64)
        D3, D4 = edgeworth(mu, C, F["k3"][s_].astype(np.float64), F["k4"][s_].astype(np.float64), D21, K22, K31)
        res[h] = dict(g=g, e3=qdiag(W[t], D3), e4=qdiag(W[t], D4), vt=F["var"][t].astype(np.float64))
    ref = res["full"]["vt"]
    A, B = res["h0"], res["h1"]

    def nf(fa, fb):
        return np.sqrt(max(np.mean((fa / ref) * (fb / ref)), 0.0))

    def share(fa, fb, ya, yb):
        """noise-corrected share of y's energy explained by projection on f (cross halves)"""
        num = np.mean(fa * yb / ref ** 2) * np.mean(fb * ya / ref ** 2)
        den = np.mean(fa * fb / ref ** 2) * np.mean(ya * yb / ref ** 2)
        return num / den if den > 0 else np.nan

    print(f" layer {t}:", flush=True)
    print(f"   g (true defect)   mean {np.mean(res['full']['g'] / ref):+.2e}  nf {nf(A['g'], B['g']):.2e}", flush=True)
    print(f"   E3 (order 3)      mean {np.mean(res['full']['e3'] / ref):+.2e}  nf {nf(A['e3'], B['e3']):.2e}", flush=True)
    print(f"   E4 (order 4)      mean {np.mean(res['full']['e4'] / ref):+.2e}  nf {nf(A['e4'], B['e4']):.2e}", flush=True)
    rA, rB = A["g"] - A["e3"] - A["e4"], B["g"] - B["e3"] - B["e4"]
    r3A, r3B = A["g"] - A["e3"], B["g"] - B["e3"]
    print(f"   g - E3            mean {np.mean((res['full']['g'] - res['full']['e3']) / ref):+.2e}  nf {nf(r3A, r3B):.2e}", flush=True)
    print(f"   g - E3 - E4       mean {np.mean((res['full']['g'] - res['full']['e3'] - res['full']['e4']) / ref):+.2e}  nf {nf(rA, rB):.2e}", flush=True)
    print(f"   x (chain)         mean {np.mean(x / ref):+.2e}  rms {np.sqrt(np.mean((x / ref) ** 2)):.2e}", flush=True)
    xa, xb = x - A["g"], x - B["g"]
    print(f"   x - g (injected)  mean {np.mean((x - res['full']['g']) / ref):+.2e}  nf {nf(xa, xb):.2e}", flush=True)
    xe_a, xe_b = x - A["e3"] - A["e4"], x - B["e3"] - B["e4"]
    print(f"   x - E3 - E4       mean {np.mean((x - res['full']['e3'] - res['full']['e4']) / ref):+.2e}  nf {nf(xe_a, xe_b):.2e}", flush=True)
    print(f"   share of (x - g) explained by -(g - E3 - E4): {share(-rA, -rB, xa, xb):+.2f};  by (x - E3 - E4): "
          f"{share(xe_a, xe_b, xa, xb):+.2f}", flush=True)
    for k in ("g", "e3", "e4"):
        for h, tag in (("full", "F"), ("h0", "A"), ("h1", "B")):
            save[f"{k}{tag}_{t}"] = res[h][k]
    save[f"x_{t}"] = x; save[f"vt_{t}"] = ref
np.savez(f"var_ladder_{net}.npz", **save)
