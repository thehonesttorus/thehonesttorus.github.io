# The chain's per-layer statistics against Monte Carlo truth (note XXXI).   python oracle_cmp.py NET PREFIX [onestep [TRUTH]]
# With "onestep": reads chain_off{NET}_o.npz (a dump run with every oracle on) and compares the chain's own value of each
# statistic at layer l, computed from Monte Carlo truth at layer l - 1 before the oracle replaced it ("<name>own").
# TRUTH (default full) names the Monte Carlo file used as truth: run the oracles on h0 and judge against h1 so that the
# inputs' and the target's sampling noise are independent.
# For every statistic the oracle hooks can replace: relative error ||chain - MC|| / ||MC||, correlation, and the Monte
# Carlo noise of the full-sample estimate (half the half-sample difference over ||MC||). Off-diagonal entries only for
# the matrices. Chain names: var, D3, D21, g4row (kappa_4 diagonal), wk4m (kappa(z_i,z_i,z_j,z_j)), wk431
# (wk431[a, c] = kappa(z_a, z_c, z_c, z_c), compared with the transposed MC (3,1) slice) and C_off.
import sys, numpy as np
net, pre = int(sys.argv[1]), sys.argv[2]
one = len(sys.argv) > 3 and sys.argv[3] == "onestep"; sfx = "own" if one else ""
tt = sys.argv[4] if len(sys.argv) > 4 else "full"
ch = np.load(f"chain_off{net}{'_o' if one else ''}.npz"); F = np.load(f"{pre}_{tt}.npz"); H0 = np.load(f"{pre}_h0.npz"); H1 = np.load(f"{pre}_h1.npz")
print(f"net {net}: truth {tt}, {int(F['n'])} samples" + ("; one-step: chain's own value from true inputs at the previous layer" if one else ""))
off = ~np.eye(1024, dtype=bool)
PAIRS = [("var", "var", False, False), ("D3", "k3", False, False), ("D21", "D21", True, False), ("g4row", "k4", False, False),
         ("wk4m", "K22", True, False), ("wk431", "K31", True, True), ("C_off", "cov", True, False)]
for l in range(1, 16):
    parts = []
    for cn, mn, isoff, tr in PAIRS:
        if f"{cn}{sfx}_{l}" not in ch.files or mn not in F.files:
            continue
        c = ch[f"{cn}{sfx}_{l}"].astype(np.float64)
        t, h0, h1 = (X[mn][l].astype(np.float64) for X in (F, H0, H1))
        if tr:
            t, h0, h1 = t.T, h0.T, h1.T
        if isoff:
            c, t, h0, h1 = c[off], t[off], h0[off], h1[off]
        nt = np.linalg.norm(t)
        parts.append(f"{cn} rel {np.linalg.norm(c - t) / nt:.3f} corr {np.corrcoef(c, t)[0, 1]:.4f} "
                     f"(MC {np.linalg.norm(h0 - h1) / 2 / nt:.3f})")
    if parts:
        print(f"layer {l:2d}: " + " | ".join(parts), flush=True)
