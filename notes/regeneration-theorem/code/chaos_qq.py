# Note XLIX section 9: the cross-age second-chaos pair term of the (3,1) slice, from the chaos state.
#   kappa(z_i, z_i, z_i, z_j) >= 6 sum_(s,t) J2_s J2_t V_is V_it U_is G_st U_jt + 6 sum_(s,t) J2_s J2_t V_is V_jt U_is G_st U_it
# over source pairs born at DIFFERENT layers (same-layer pairs are the ladder's C-trees). Sources: folds at every layer
# m <= L, birth direction l_s = row of L^m, jet J2_s (true marginal with kappa3), leg = column of W_(m+1) transported by
# first jets; U = L^(L+1) Lambda, G_st = l_s . l_t. Combined with the saved ladder of section 8.
#   python chaos_qq.py NET L MCDIR WDIR LADDER_NPZ
import sys, os, re, math, numpy as np
from scipy.special import ndtr
net, L, MC, WD, LAD = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4], sys.argv[5]
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "regen4_fast.py")).read()
g = {"np": np, "math": math, "re": re, "ndtr": ndtr, "SQ2PI": math.sqrt(2 * math.pi), "EMAX": 4}
f64 = lambda a: np.asarray(a, dtype=np.float64)
T = np.load(f"{MC}/mc2_off{net}_full.npz")
Wall = np.load(f"{WD}/W_off{net}.npy", mmap_mode="r")
n = Wall.shape[1]; g["n"] = n
exec(src[src.index("def he("):src.index("def r2(")], g)
mu, var, k3 = f64(T["mu"]), f64(T["var"]), f64(T["k3"])
Lm = f64(Wall[0])
blocks = []          # per birth layer m: (Lam_m (n_in x n), J2_m (n), V_m (n x n) transported to the current layer)
old_at_fold = []     # per fold layer m: the legs of every older source as seen by z^m (list over older blocks)
for m in range(L + 1):
    F = g["site_factors"](mu[m], var[m], k3[m], 0 * mu[m], True)
    Phi, J2 = F[(1, 1)], F[(1, 2)]
    Wn = f64(Wall[m + 1])
    old_at_fold.append([Vb.copy() for (_, _, Vb) in blocks])
    blocks = [(Lb, Jb, Wn @ (Phi[:, None] * Vb)) for (Lb, Jb, Vb) in blocks]
    blocks.append((Lm.T.copy(), J2, Wn.copy()))
    Lm = Wn @ (Phi[:, None] * Lm)
# U blocks at layer L+1
Ub = [Lm @ Lb for (Lb, _, _) in blocks]
Ab = [Vb * U * Jb[None, :] for (_, Jb, Vb), U in zip(blocks, Ub)]
Ptot = sum(A @ Lb.T for A, (Lb, _, _) in zip(Ab, blocks))      # n x n_in
K_cross = np.zeros((n, n)); K_same = np.zeros((n, n))
for A, (Lb, Jb, Vb), U in zip(Ab, blocks, Ub):
    Psame = A @ Lb.T
    for P, K in ((Ptot - Psame, K_cross), (Psame, K_same)):
        Pi = P @ Lb                                    # n x n: Pi_it = sum_s A_is G_st over the chosen s
        K += 6.0 * ((Pi * Vb * Jb[None, :]) @ U.T) + 6.0 * ((Pi * U * Jb[None, :]) @ Vb.T)
# secondary third chaos: the fold at (m, a) acting on the second chaos z^m_a already carries,
# kernel J2_a Sym(l_a x A_a), A_a = (1/2) sum_(s older) J2_s V^(m)_as l_s l_s^T; read by (T, L, L, L) Wick pairings
K_sec = np.zeros((n, n))
for m in range(1, L + 1):
    Ut, Jm, Vt = Ub[m], blocks[m][1], blocks[m][2]          # fold layer m: arms U~ (= U of block m), jets, legs
    q = np.zeros((n, n))
    for b, Vold in enumerate(old_at_fold[m]):
        q += (Ub[b] ** 2 * blocks[b][1][None, :]) @ Vold.T     # q_ia = sum_s V^(m)_as J2_s U_is^2
        R = (Vt * Ut * (2.0 * Jm)[None, :]) @ Vold              # R_is = sum_a V''_ia 2 J2_a U~_ia V^(m)_as
        K_sec += 3.0 * ((R * Ub[b] * blocks[b][1][None, :]) @ Ub[b].T)
    K_sec += 3.0 * ((Vt * Jm[None, :] * q) @ Ut.T) + 3.0 * ((Ut * q * Jm[None, :]) @ Vt.T)
d = np.load(LAD)
mask = d["mask"].astype(bool); w = d["w"]; t0, t1, tf, cum = d["t0"], d["t1"], d["tf"], d["cum"]
def r2v(y):
    ip = lambda A, B: float(np.sum(A * B * w * w))
    return 1.0 - ip(t0 - y, t1 - y) / ip(t0, t1), ip(y, tf) / ip(y, y)
yc, ys = K_cross[mask], K_same[mask]
print(f"=== network {net}, (3,1) slice of layer {L + 1}: second-chaos pair terms from the chaos state ({len(blocks)} birth layers)")
ysec = K_sec[mask]
for name, y in (("ladder (section 8)", cum), ("QQ cross-age alone", yc), ("QQ same-age alone", ys),
                ("secondary stars alone", ysec), ("ladder + QQ cross-age", cum + yc),
                ("ladder + QQ cross + sec. stars", cum + yc + ysec)):
    R, s = r2v(y)
    print(f"  {name:24s} R2 {100 * R:6.1f}%  scale {s:+.3f}")
ip = lambda A, B: float(np.sum(A * B * w * w))
print(f"  energy ratio |QQ cross| / |ladder| = {math.sqrt(ip(yc, yc) / ip(cum, cum)):.3f};  "
      f"corr(QQ same, ladder) = {ip(ys, cum) / math.sqrt(ip(ys, ys) * ip(cum, cum)):+.3f}")
