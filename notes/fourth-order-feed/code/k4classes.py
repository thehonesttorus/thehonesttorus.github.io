# The fourth cumulant of the pre-activation z_{l+1} = W h_l decomposed by index-multiplicity class of the post-activation
# indices, without any fourth-order tensor: per sample and output neuron a_j = W_ij X_j (X = h_l - E h_l), p_r = sum_j a_j^r,
#   class (4): p4;  (3+1): 4(p1 p3 - p4);  (2+2): 3(p2^2 - p4);  (2+1+1): 6(p1^2 p2 - 2 p1 p3 - p2^2 + 2 p4);
#   (1+1+1+1): p1^4 - 6 p1^2 p2 + 3 p2^2 + 8 p1 p3 - 6 p4   (sum = p1^4 = z^4).
# The Gaussian pairings of each class are subtracted with the sample covariance C of h_l (inclusion-exclusion over
# coincident indices), so the classes sum to kappa_4(z_i) = E z^4 - 3 (E z^2)^2 exactly.  Each class is then projected on
# the scale-mixture pattern 3 sigma_z^4 (the g_4 of note XIX).  Classes (4) + (2+2) are what the pair-slice ledger keeps.
# One Monte Carlo pass of T inputs through the official network, accumulating at the requested source layers.
import numpy as np, sys, time
net = int(sys.argv[1]) if len(sys.argv) > 1 else 0
layers = [int(v) for v in (sys.argv[2] if len(sys.argv) > 2 else "7,10,14").split(",")]   # source layers l (target z_{l+1})
T = int(sys.argv[3]) if len(sys.argv) > 3 else 400000; B = 2000
Wc = np.load(f"../official/W_off{net}.npy").astype(np.float64); L, n = Wc.shape[0], Wc.shape[1]
D = np.load(f"mc_cum_off{net}.npz"); m_h = {l: D["s1h"][l].astype(np.float64) for l in layers}
acc = {l: dict(P4=np.zeros(n), P13=np.zeros(n), P22=np.zeros(n), P112=np.zeros(n), P1111=np.zeros(n), P11=np.zeros(n), C=np.zeros((n, n))) for l in layers}
Wp = {l: [Wc[l + 1] ** r for r in (1, 2, 3, 4)] for l in layers}
rng = np.random.default_rng(1234 + net); t0 = time.time(); done = 0
while done < T:
    nb = min(B, T - done); h = rng.standard_normal((nb, n))
    for l in range(max(layers) + 1):
        y = h @ Wc[l].T; h = np.maximum(y, 0.0)
        if l in layers:
            X = h - m_h[l]; a = acc[l]; a["C"] += X.T @ X
            P1 = X @ Wp[l][0].T; P2 = (X * X) @ Wp[l][1].T; P3 = (X ** 3) @ Wp[l][2].T; P4 = (X ** 4) @ Wp[l][3].T
            a["P4"] += P4.sum(0); a["P13"] += (P1 * P3).sum(0); a["P22"] += (P2 * P2).sum(0); a["P112"] += (P1 * P1 * P2).sum(0); a["P1111"] += (P1 ** 4).sum(0); a["P11"] += (P1 * P1).sum(0)
    done += nb
    if done % 40000 == 0: print(f"  {done}/{T} samples, {time.time()-t0:.0f}s", flush=True)
def fit(x, y): return (x @ y) / (x @ x)
print(f"net {net}, T = {T}: kappa_4 of z_(l+1) by index-multiplicity class of h_l, projected on 3 sigma_z^4 (g_4 by class); 'pair' = (4) + (2+2), the ledger's retained part")
print("  l+1 | g4 total | (4) | (3+1) | (2+2) | (2+1+1) | (1+1+1+1) | pair | pair/total | sigma check | R^2 of total fit")
out = {}
for l in layers:
    a = acc[l]; W = Wc[l + 1]; C = a["C"] / T
    E4, E13, E22, E112, E1111, E11 = (a[k] / T for k in ("P4", "P13", "P22", "P112", "P1111", "P11"))
    R4 = E4; R31 = 4 * (E13 - E4); R22 = 3 * (E22 - E4); R211 = 6 * (E112 - 2 * E13 - E22 + 2 * E4); R1111 = E1111 - 6 * E112 + 3 * E22 + 8 * E13 - 6 * E4
    c = np.diag(C); W2 = W * W; W3 = W2 * W; W4 = W2 * W2
    WC = W @ C; sig2 = np.einsum("ij,ij->i", WC, W)           # full variance of z_i (from the same samples' covariance)
    Dg = W2 @ c                                                # sum_j W_ij^2 C_jj
    Q = np.einsum("ij,ij->i", W2 @ (C * C), W2)                # sum_jk W_ij^2 W_ik^2 C_jk^2 (all j,k)
    G4 = 3 * W4 @ (c * c)
    G31 = 12 * np.einsum("ij,ij->i", W3 * c[None, :], WC - W * c[None, :])
    G22 = 3 * (Dg * Dg - W4 @ (c * c)) + 6 * (Q - W4 @ (c * c))
    # (2+1+1): sum over distinct (j,k,l) of W_ij^2 W_ik W_il [C_jj C_kl + 2 C_jk C_jl], times 6
    t1_all = Dg * sig2; t1_kj = np.einsum("ij,ij->i", W2 * W * c[None, :], WC); t1_kl = Dg * Dg; t1_jjj = W4 @ (c * c)
    t1 = t1_all - 2 * t1_kj - t1_kl + 2 * t1_jjj
    t2_all = np.einsum("ij,ij->i", W2, WC * WC); t2_kj = np.einsum("ij,ij->i", W2 * W * c[None, :], WC); t2_kl = Q; t2_jjj = W4 @ (c * c)
    t2 = t2_all - 2 * t2_kj - t2_kl + 2 * t2_jjj
    G211 = 6 * (t1 + 2 * t2)
    G1111 = 3 * sig2 * sig2 - G4 - G31 - G22 - G211
    K = dict(k4=R4 - G4, k31=R31 - G31, k22=R22 - G22, k211=R211 - G211, k1111=R1111 - G1111)
    tot = sum(K.values()); x4 = 3 * E11 * E11; gt = fit(x4, tot); r2 = 1 - np.var(tot - gt * x4) / np.var(tot)
    g = {k: fit(x4, v) for k, v in K.items()}; pair = g["k4"] + g["k22"]
    print(f"  {l+1:3d} | {gt:+.5f} | {g['k4']:+.5f} | {g['k31']:+.5f} | {g['k22']:+.5f} | {g['k211']:+.5f} | {g['k1111']:+.5f} | {pair:+.5f} | {pair/gt:.3f} | {np.mean(sig2/E11):.4f} | {r2:.2f}", flush=True)
    out[l] = dict(g=g, gt=gt, sig2=sig2, E11=E11, K=K)
    # the mixture's prediction for the dropped classes: (2+1+1) ~ 6 g sigma_diag^2 sigma_off^2, (1111) ~ 3 g sigma_off^4 at leading order
    off = sig2 - Dg; print(f"        off-diagonal variance share mean {np.mean(off/sig2):+.4f}; mixture prediction for dropped classes at g = {gt:.4f}: (3+1) ~ {fit(x4, 12*gt*np.einsum('ij,ij->i', W3*c[None,:], WC - W*c[None,:])):+.5f}, (2+1+1) ~ {fit(x4, 6*gt*Dg*off):+.5f}, (1111) ~ {fit(x4, 3*gt*off*off):+.5f}")
np.savez(f"k4classes_off{net}.npz", **{f"l{l}_{k}": v for l in layers for k, v in out[l]["g"].items()}, **{f"l{l}_gt": out[l]["gt"] for l in layers})
