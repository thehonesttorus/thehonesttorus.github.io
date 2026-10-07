# The one-loop (Gaussian reference) fourth cumulant of z_{l+1} = W relu(y), y ~ N(mu_l, S_l) with the TRUE marginal state
# of layer l, decomposed by index-multiplicity class exactly as k4classes.py does for the true law.  The difference
# true - Gaussian is the transported content at first order; the mixture prediction for its dropped classes is
# (2+1+1) ~ 6 g4(l) sigma_diag^2 sigma_off^2 and (1111) ~ 3 g4(l) sigma_off^4 with g4(l) the measured gain of layer l.
import numpy as np, sys, time
from scipy.special import ndtr
net = int(sys.argv[1]) if len(sys.argv) > 1 else 0
layers = [int(v) for v in (sys.argv[2] if len(sys.argv) > 2 else ",".join(str(i) for i in range(1, 15))).split(",")]
T = int(sys.argv[3]) if len(sys.argv) > 3 else 300000; B = 2000
Wc = np.load(f"../official/W_off{net}.npy").astype(np.float64); L, n = Wc.shape[0], Wc.shape[1]
D = np.load(f"mc_cum_off{net}.npz"); Z = np.load(f"gac3_off{net}.npz"); gm4 = Z["gm4"]
def fit(x, y): return (x @ y) / (x @ x)
print(f"net {net}, T = {T}: one-loop kappa_4 of z_(l+1) by class (Gaussian y_l with the true marginal state), projected on 3 sigma_z^4")
print("  l+1 | phi4 total | (4) | (3+1) | (2+2) | (2+1+1) | (1+1+1+1) | pair | off share (sigma^4-weighted) | mixture (2+1+1) at g4(l) | mixture (1111)")
rng = np.random.default_rng(777 + net); t0 = time.time()
for l in layers:
    mu = D["s1y"][l].astype(np.float64); S = D["S2y"][l].astype(np.float64) - np.outer(mu, mu); sig = np.sqrt(np.diag(S)); a = mu / sig
    mG = sig * (a * ndtr(a) + np.exp(-a * a / 2) / np.sqrt(2 * np.pi))
    ev, U = np.linalg.eigh(S); ev = np.maximum(ev, 1e-10 * ev.max()); Lc = U * np.sqrt(ev)
    W = Wc[l + 1]; Wp = [W ** r for r in (1, 2, 3, 4)]
    acc = dict(P4=np.zeros(n), P13=np.zeros(n), P22=np.zeros(n), P112=np.zeros(n), P1111=np.zeros(n), P11=np.zeros(n), C=np.zeros((n, n)), m=np.zeros(n))
    done = 0
    while done < T:
        nb = min(B, T - done); y = mu + rng.standard_normal((nb, n)) @ Lc.T; h = np.maximum(y, 0.0); X = h - mG
        acc["m"] += X.sum(0); acc["C"] += X.T @ X
        P1 = X @ Wp[0].T; P2 = (X * X) @ Wp[1].T; P3 = (X ** 3) @ Wp[2].T; P4 = (X ** 4) @ Wp[3].T
        acc["P4"] += P4.sum(0); acc["P13"] += (P1 * P3).sum(0); acc["P22"] += (P2 * P2).sum(0); acc["P112"] += (P1 * P1 * P2).sum(0); acc["P1111"] += (P1 ** 4).sum(0); acc["P11"] += (P1 * P1).sum(0)
        done += nb
    C = acc["C"] / T - np.outer(acc["m"] / T, acc["m"] / T)
    E4, E13, E22, E112, E1111, E11 = (acc[k] / T for k in ("P4", "P13", "P22", "P112", "P1111", "P11"))
    R4 = E4; R31 = 4 * (E13 - E4); R22 = 3 * (E22 - E4); R211 = 6 * (E112 - 2 * E13 - E22 + 2 * E4); R1111 = E1111 - 6 * E112 + 3 * E22 + 8 * E13 - 6 * E4
    c = np.diag(C); W2 = W * W; W3 = W2 * W; W4 = W2 * W2; WC = W @ C; sig2 = np.einsum("ij,ij->i", WC, W); Dg = W2 @ c
    Q = np.einsum("ij,ij->i", W2 @ (C * C), W2)
    G4 = 3 * W4 @ (c * c); G31 = 12 * np.einsum("ij,ij->i", W3 * c[None, :], WC - W * c[None, :]); G22 = 3 * (Dg * Dg - W4 @ (c * c)) + 6 * (Q - W4 @ (c * c))
    t1 = Dg * sig2 - 2 * np.einsum("ij,ij->i", W2 * W * c[None, :], WC) - Dg * Dg + 2 * W4 @ (c * c)
    t2 = np.einsum("ij,ij->i", W2, WC * WC) - 2 * np.einsum("ij,ij->i", W2 * W * c[None, :], WC) - Q + 2 * W4 @ (c * c)
    G211 = 6 * (t1 + 2 * t2); G1111 = 3 * sig2 * sig2 - G4 - G31 - G22 - G211
    K = dict(k4=R4 - G4, k31=R31 - G31, k22=R22 - G22, k211=R211 - G211, k1111=R1111 - G1111)
    x4 = 3 * E11 * E11; g = {k: fit(x4, v) for k, v in K.items()}; gt = sum(g.values()); off = sig2 - Dg
    offw = fit(x4, 3 * sig2 * off) / 1.0   # sigma^4-weighted mean of off/sig2 (projection of 3 sig2 off on 3 sig2^2)
    print(f"  {l+1:3d} | {gt:+.5f} | {g['k4']:+.5f} | {g['k31']:+.5f} | {g['k22']:+.5f} | {g['k211']:+.5f} | {g['k1111']:+.5f} | {g['k4']+g['k22']:+.5f} | {offw:+.4f} | {fit(x4, 6*gm4[l]*Dg*off):+.5f} | {fit(x4, 3*gm4[l]*off*off):+.5f}   [{time.time()-t0:.0f}s]", flush=True)
    np.savez(f"k4gauss_off{net}_l{l}.npz", **g, gt=gt, offw=offw, m211=fit(x4, 6*gm4[l]*Dg*off), m1111=fit(x4, 3*gm4[l]*off*off))
