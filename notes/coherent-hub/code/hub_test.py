# Note XLV section 3: one-step regeneration test of the coherent doubly-occupied hub terms.
#   python -I hub_test.py HELPERS.py LOCDIR MC2DIR CD2DIR NET
# HELPERS.py = phaseM/var_ladder.py (only its definitions are used: ccoef, series, relu_var, qdiag).
# From the true layer s-1 state (mu, C) and W_s:
#   P = W diag(Phi) C_off = Cov(z', z) (first jet),  G = 2 [Phi (1 - Phi) - alpha phi Phi - phi^2],  Q = 2 m (1 - Phi)
#   H31_ij = 3 sum_a W_ia^2 G_a P_ia P_ja,  H22_ij = sum_a G_a (W_ia^2 P_ja^2 + W_ja^2 P_ia^2),  H4_i = 6 sum_a W_ia^2 G_a P_ia^2
#   HD21_ij = sum_a W_ia^2 Q_a P_ja,  HD3_i = 3 sum_a W_ia^2 Q_a P_ia
#   annealed forms: W_ia^2 -> mean(W^2).
# Compared with the true layer-s slices (K31[i,j] = kappa(i,i,i,j), K22, k4, D21[i,j] = kappa(i,i,j), k3) and the chain's
# own (chaindump2), on neurons with true alpha > -2.5, entrywise and in the slice's readout into layer s+1's variance.
import sys, numpy as np
from scipy.special import ndtr

hp, loc, mcd, cdd, net = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5])
_src = open(hp).read().split("loc, mcd, net = sys.argv")[0] + open(hp).read().split("layers = [int(v) for v in sys.argv[4].split(\",\")]")[1].split("cd = np.load")[0]
exec(_src)
SQ2PI = np.sqrt(2.0 * np.pi)
f64 = lambda a: np.asarray(a, dtype=np.float64)
sym = lambda A: 0.5 * (A + A.T)
W = np.load(f"{loc}/W_off{net}.npy").astype(np.float64)
c2 = np.load(f"{cdd}/chaindump2_{net}.npz")
T = {h: np.load(f"{mcd}/mc2_off{net}_{h}.npz") for h in ("full", "h0", "h1")}
n = W.shape[1]
off = ~np.eye(n, dtype=bool)


def zd(A):
    A = A.copy(); np.fill_diagonal(A, 0.0); return A


def hubs(mu, C, Ws):
    v = np.diag(C).copy(); S = np.sqrt(v); al = mu / S
    Phi = ndtr(al); phi = np.exp(-0.5 * al * al) / SQ2PI; m = S * (al * Phi + phi)
    G = 2.0 * (Phi * (1.0 - Phi) - al * phi * Phi - phi * phi); Q = 2.0 * m * (1.0 - Phi)
    P = Ws @ (Phi[:, None] * zd(C))
    WW = Ws * Ws; mw = float(WW.mean())
    u = (P * P) @ G
    H = {}
    H["K31q"] = zd(3.0 * ((WW * P * G[None, :]) @ P.T))
    H["K31a"] = zd(3.0 * mw * ((P * G[None, :]) @ P.T))
    A = (WW * G[None, :]) @ (P * P).T
    H["K22q"] = zd(A + A.T)
    H["K22a"] = zd(mw * (u[:, None] + u[None, :]))
    H["k4q"] = 6.0 * ((WW * P * P) @ G)
    H["k4a"] = 6.0 * mw * u
    H["D21q"] = zd((WW * Q[None, :]) @ P.T)
    H["D21a"] = zd(mw * (np.ones((n, 1)) @ (P @ Q)[None, :]))
    H["k3q"] = 3.0 * ((WW * P) @ Q)
    H["k3a"] = 3.0 * mw * (P @ Q)
    return H


def readout_map(mu, C):
    """linear maps from each slice to its contribution to Cov(relu(z)) at the layer (exact bivariate Edgeworth)."""
    v = np.diag(C).copy(); s = np.sqrt(v); al = mu / s
    R = C / np.outer(s, s); np.fill_diagonal(R, 0.0)
    c = ccoef(al, M + 5); sa = s[:, None]; sb = s[None, :]
    S31 = series(c, c, R, 3, 1, 0) / 6.0 * sa ** -2
    S22 = 0.25 * series(c, c, R, 2, 2, 0) / (sa * sb)
    S21 = 0.5 * series(c, c, R, 2, 1, 0) / sa
    S40 = series(c, c, R, 4, 0, 1) / 24.0 * sa ** -3 * sb
    S30 = series(c, c, R, 3, 0, 1) / 6.0 * sa ** -2 * sb
    _, m1 = relu_var(mu, v)
    d4 = (2.0 * c[3] / s ** 2 - 2.0 * m1 * c[4] / s ** 3) / 24.0
    d3 = (2.0 * c[2] / s - 2.0 * m1 * c[3] / s ** 2) / 6.0

    def K31(X):
        Y = X * S31; Y = Y + Y.T; np.fill_diagonal(Y, 0.0); return Y

    def K22(X):
        Y = X * S22; np.fill_diagonal(Y, 0.0); return Y

    def D21(X):
        Y = X * S21; Y = Y + Y.T; np.fill_diagonal(Y, 0.0); return Y

    def k4(x):
        Y = x[:, None] * S40; Y = Y + Y.T; np.fill_diagonal(Y, x * d4); return Y

    def k3(x):
        Y = x[:, None] * S30; Y = Y + Y.T; np.fill_diagonal(Y, x * d3); return Y
    return dict(K31=K31, K22=K22, D21=D21, k4=k4, k3=k3)


def truth_slice(h, key, s):
    F = T[h]
    if key == "K31":
        return zd(f64(F["K31"][s]))
    if key == "K22":
        return zd(sym(f64(F["K22"][s])))
    if key == "D21":
        return zd(f64(F["D21"][s]))
    return f64(F[key][s])


def chain_slice(key, s):
    if key == "K31":
        return zd(f64(c2[f"wk431_{s}"]).T)
    if key == "K22":
        return zd(sym(f64(c2[f"wk4m_{s}"])))
    if key == "D21":
        return zd(f64(c2[f"D21_{s}"]))
    if key == "k4":
        return f64(c2[f"g4row_{s}"])
    return f64(c2[f"D3_{s}"])


print(f"=== network {net}: coherent doubly-occupied hub, one-step regeneration from the true layer s-1 state", flush=True)
print("per slice: fit truth ~ c_p package + c_h H (c_h, R2 package alone -> joint);  chain error ~ c H: coefficient and "
      "share of the chain-error energy, entrywise | in the next-layer variance readout;  q = quenched, a = annealed", flush=True)
for s in (4, 6, 8, 10, 12, 13, 14):
    Ff = T["full"]
    mu0 = f64(Ff["mu"][s - 1]); C0 = sym(f64(Ff["cov"][s - 1]))
    H = hubs(mu0, C0, W[s])
    mu1 = f64(Ff["mu"][s]); C1 = sym(f64(Ff["cov"][s])); v1 = np.diag(C1).copy(); C1o = zd(C1)
    act = (mu1 / np.sqrt(v1)) > -2.5
    mk = off & act[:, None] & act[None, :]
    mu2 = f64(Ff["mu"][s + 1]); v2 = f64(Ff["var"][s + 1]); act2 = (mu2 / np.sqrt(v2)) > -2.5
    forms = dict(K31=12.0 * v1[:, None] * C1o, K22=zd(4.0 * np.outer(v1, v1) + 8.0 * C1o * C1o), k4=12.0 * v1 * v1,
                 D21=zd(2.0 * (2.0 * mu1[:, None] * C1 + mu1[None, :] * v1[:, None])), k3=6.0 * mu1 * v1)
    RM = readout_map(mu1, C1)
    Wn = W[s + 1]
    line = [f"  s={s:2d}"]
    for key in ("K31", "K22", "k4", "D21", "k3"):
        mat = key in ("K31", "K22", "D21")
        msk = mk if mat else act
        ip = lambda A, B: float(np.sum((A * B)[msk]))
        Tf, T0, T1 = (truth_slice(h, key, s) for h in ("full", "h0", "h1"))
        f = forms[key]
        res = []
        for form in ("q", "a"):
            Hx = H[key + form]
            # joint fit on the full sample, residual energy across halves
            A2 = np.array([[ip(f, f), ip(f, Hx)], [ip(f, Hx), ip(Hx, Hx)]]); b2 = np.array([ip(f, Tf), ip(Hx, Tf)])
            cp, ch = np.linalg.solve(A2, b2)
            cp0 = ip(f, Tf) / ip(f, f)
            ET = ip(T0, T1)
            r2p = 1.0 - ip(T0 - cp0 * f, T1 - cp0 * f) / ET
            r2j = 1.0 - ip(T0 - cp * f - ch * Hx, T1 - cp * f - ch * Hx) / ET
            cs = ""
            if s <= 14:
                Cc = chain_slice(key, s)
                e0, e1, ef = T0 - Cc, T1 - Cc, Tf - Cc
                EE = ip(e0, e1)
                c = ip(ef, Hx) / ip(Hx, Hx)
                sh = 1.0 - ip(e0 - c * Hx, e1 - c * Hx) / EE
                # readout metric: contribution to the next layer's variance (active next-layer neurons)
                ro = lambda X: qdiag(Wn, RM[key](np.where(msk, X, 0.0)))[act2]
                r0, r1, rf, rh = ro(e0), ro(e1), ro(ef), ro(Hx)
                cr = float(rf @ rh / (rh @ rh))
                shr = 1.0 - float((r0 - cr * rh) @ (r1 - cr * rh)) / float(r0 @ r1)
                cs = f" err c {c:+.2f} sh {100 * sh:5.1f}% | ro c {cr:+.2f} sh {100 * shr:5.1f}%"
            res.append(f"{form}: c_h {ch:+.2f} R2 {100 * r2p:5.1f}->{100 * r2j:5.1f}{cs}")
        line.append(f"    {key:3s} " + "  ||  ".join(res))
    print("\n".join(line), flush=True)
