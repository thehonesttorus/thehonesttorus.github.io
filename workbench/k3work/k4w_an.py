# The (2,2) pre-activation slice: chain's wk4m against Monte Carlo truth, and candidate constructions.
import sys, pickle, numpy as np
net = int(sys.argv[1]); srcs = [int(v) for v in sys.argv[2].split(",")]
Z = np.load(f"k4mc2_off{net}.npz"); T = float(Z["T"]); n = 1024
Wc = np.load(f"../official/W_off{net}.npy").astype(np.float64)
D = {d["layer"]: d for d in pickle.load(open(f"v29k4_off{net}.pkl", "rb"))["dumps"]}
iu = np.triu_indices(n, 1)
def cum(kind, l, hf):
    g = (lambda k: Z[f"h{hf}_{kind}{l}_{k}"] / T) if hf is not None else (lambda k: 0.5 * (Z[f"h0_{kind}{l}_{k}"] + Z[f"h1_{kind}{l}_{k}"]) / T)
    x1 = g("s1"); C = g("C") - np.outer(x1, x1); c = np.diag(C).copy()
    K22 = g("M22") - np.outer(c, c) - 2 * C * C; np.fill_diagonal(K22, 0.0)
    return C, c, K22, g("s4") - 3 * c * c
def stat(est, tru, noise):
    a, b = est[iu], tru[iu]; coef = float(a @ b / (a @ a))
    return f"corr {np.corrcoef(a, b)[0, 1]:+.3f} relerr {np.linalg.norm(a - b) / np.linalg.norm(b):.3f} LScoef {coef:.3f} (truth noise rel {np.linalg.norm(noise[iu]) / np.linalg.norm(b):.3f})"
for s in srcs:
    t = s + 1; W = Wc[t]; W2 = W * W
    Cy, cy, K22y, _ = cum("y", t, None); _, _, K22y0, _ = cum("y", t, 0); _, _, K22y1, _ = cum("y", t, 1)
    noise = 0.5 * (K22y0 - K22y1)
    Cx, cx, K22x, K4x = cum("h", s, None)
    wk = D[t]["wk4m"].astype(np.float64); g4 = D[t]["g4row"].astype(np.float64); var = D[t]["var"].astype(np.float64)
    Coff = D[t]["C_off"].astype(np.float64)
    print(f"\n source {s} -> target {t}: truth K22(y) offdiag rms {np.sqrt(np.mean(K22y[iu]**2)):.3e}, mean {K22y[iu].mean():+.3e}")
    print("  chain wk4m (g4_i + g4_j)/6        " + stat(wk, K22y, noise))
    # pair class from the TRUE post-activation slices, exact in the weights (diagonal class + (2,2) class, no cross term)
    P1 = W2 @ K22x @ W2.T + W2 @ (K4x[:, None] * W2.T)
    np.fill_diagonal(P1, 0.0)
    print("  pair class, true K22(x), exact W   " + stat(P1, K22y, noise))
    # the same from the chain's own post-activation slices
    K22c = D[s]["K22"].astype(np.float64).copy(); np.fill_diagonal(K22c, 0.0); K22c = 0.5 * (K22c + K22c.T); K4c = D[s]["K4v"].astype(np.float64)
    P2 = W2 @ K22c @ W2.T + W2 @ (K4c[:, None] * W2.T); np.fill_diagonal(P2, 0.0)
    print("  pair class, chain K22(x), exact W  " + stat(P2, K22y, noise))
    lam, V = np.linalg.eigh(K22c); o = np.argsort(-np.abs(lam))
    for k in (4, 16):
        Vk, lk = V[:, o[:k]], lam[o[:k]]; WV = W2 @ Vk
        P3 = (WV * lk) @ WV.T + W2 @ (K4c[:, None] * W2.T); np.fill_diagonal(P3, 0.0)
        print(f"  pair class, chain K22 rank {k:2d}     " + stat(P3, K22y, noise))
    # scale-mixture product form g (var_i var_j + 2 C_ij^2), g per neuron from the chain's kappa4 diagonal
    gi = np.maximum(g4, 0) / (3 * var * var); gs = np.sqrt(np.outer(gi, gi))
    P4 = gs * (np.outer(var, var) + 2 * Coff * Coff); np.fill_diagonal(P4, 0.0)
    print("  scale mixture g(v_i v_j + 2 C^2)   " + stat(P4, K22y, noise))
    P5 = (np.add.outer(g4, g4) / 6.0) * 0 + gs * np.outer(var, var); np.fill_diagonal(P5, 0.0)
    print("  product shape g v_i v_j            " + stat(P5, K22y, noise))
    # residual structure of the chain's slice
    R = K22y - wk; np.fill_diagonal(R, 0.0)
    print("  chain residual vs 2 g C^2: " + stat(2 * gs * Coff * Coff + 0 * R, R, noise) + " | vs pair class - wk: " + stat(P2 - wk, R, noise))
