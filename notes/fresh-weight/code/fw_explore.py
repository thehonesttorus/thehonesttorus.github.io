# Note XLIV section 7 (exploratory, not pre-registered): what the collective covariance error is made of.
#   python -I fw_explore.py LOCDIR MC2DIR CD2DIR NET LADDER_PY
# Per layer s: top-k eigenbasis U of the true pre-activation covariance; B = U^T E U for the chain's error E.
#   (a) c_1, c_2, c_4, c_8 and the split of the top-8 block into eigenvalue (diagonal) and rotation energy;
#   (b) relative eigenvalue errors B_jj / lambda_j, j = 1..8 (halves averaged);
#   (c) E = P + I: P = W_s [KG(chain state s-1) - KG(true state s-1)] W_s^T is the error propagated through the
#       Gaussian closure, I = x - g the injection (chain's non-Gaussian correction minus the true one); noise-free
#       energies and top-32 shares of each.
import sys, numpy as np
_argv = list(sys.argv)
sys.argv = [_argv[5], ".", ".", "0", "1"]
exec(open(_argv[5]).read().split("cd = np.load")[0])
sys.argv = _argv
loc, mcd, cdd, net = _argv[1], _argv[2], _argv[3], int(_argv[4])
W = np.load(f"{loc}/W_off{net}.npy")
c2 = np.load(f"{cdd}/chaindump2_{net}.npz")
F = {h: np.load(f"{mcd}/mc2_off{net}_{h}.npz") for h in ("full", "h0", "h1")}
COV = {h: F[h]["cov"] for h in F}; MU = {h: F[h]["mu"] for h in F}
out = c2["out"]
n = W.shape[1]
off = ~np.eye(n, dtype=bool)
sym = lambda A: 0.5 * (A + A.T)


def chainC(s):
    C = sym(c2[f"C_off_{s}"].astype(np.float64)); np.fill_diagonal(C, c2[f"var_{s}"].astype(np.float64)); return C


trueC = lambda h, s: sym(COV[h][s].astype(np.float64))
ip = lambda A, B: float(np.sum(A[off] * B[off]))


def perp(E, U):
    EU = E @ U; UE = U.T @ E
    return E - U @ UE - EU @ U.T + U @ (U.T @ EU) @ U.T


def ck(A0, A1, U):
    return 1.0 - ip(A0, perp(A1, U)) / ip(A0, A1)


print(f"=== network {net}", flush=True)
for s in range(8, 15):
    Ws = W[s].astype(np.float64)
    lam, V = np.linalg.eigh(trueC("full", s)); lam, V = lam[::-1], V[:, ::-1]
    Cc = chainC(s)
    E = {h: Cc - trueC(h, s) for h in ("h0", "h1")}
    for h in E:
        np.fill_diagonal(E[h], 0.0)
    cs = [ck(E["h0"], E["h1"], V[:, :k]) for k in (1, 2, 4, 8)]
    U8 = V[:, :8]
    B = {h: U8.T @ E[h] @ U8 for h in E}
    dg = float(np.sum(np.diag(B["h0"]) * np.diag(B["h1"]))); tot = float(np.sum(B["h0"] * B["h1"]))
    rel = 0.5 * (np.diag(B["h0"]) + np.diag(B["h1"])) / lam[:8]
    Kc = KG(W[s - 1].astype(np.float64) @ out[s - 2].astype(np.float64), chainC(s - 1))   # chain state at s-1
    x = Cc - Ws @ Kc @ Ws.T
    P, I = {}, {}
    for h in ("h0", "h1"):
        Kt = KG(MU[h][s - 1].astype(np.float64), trueC(h, s - 1))
        g = trueC(h, s) - Ws @ Kt @ Ws.T
        P[h] = Ws @ (Kc - Kt) @ Ws.T; I[h] = x - g
        np.fill_diagonal(P[h], 0.0); np.fill_diagonal(I[h], 0.0)
    U32 = V[:, :32]
    eP, eI, eE = ip(P["h0"], P["h1"]), ip(I["h0"], I["h1"]), ip(E["h0"], E["h1"])
    cross = 0.5 * (ip(P["h0"], I["h1"]) + ip(I["h0"], P["h1"]))
    print(f"  {s:2d}  c_1/2/4/8 {cs[0]:.3f} {cs[1]:.3f} {cs[2]:.3f} {cs[3]:.3f}  top-8 block eigenvalue share {dg / tot:.3f}"
          f"  lambda_1..4/mean {lam[0] / lam.mean():.0f} {lam[1] / lam.mean():.0f} {lam[2] / lam.mean():.0f} "
          f"{lam[3] / lam.mean():.0f}", flush=True)
    print(f"      rel eigenvalue error dlam/lam j=1..8: " + " ".join(f"{r:+.1e}" for r in rel), flush=True)
    print(f"      energy split of E: propagated {eP / eE:.3f}  injected {eI / eE:.3f}  2 x cross {2 * cross / eE:+.3f};"
          f"  top-32 share: propagated {ck(P['h0'], P['h1'], U32):.3f}  injected {ck(I['h0'], I['h1'], U32):.3f}",
          flush=True)
