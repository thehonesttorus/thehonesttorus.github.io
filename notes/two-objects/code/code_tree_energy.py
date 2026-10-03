# Where does an activation's variance sit along the activation-code tree?
# Cells = sign patterns (activation codes) of k selected neurons of a given layer, chosen from the weights:
#   layer 1: neurons ranked by |<r_j, w_i>|/|w_i| (r_j = linearised row of output j; weights only)
#   layer L-1: neurons ranked by |W_L[j, i]| (the output's own incoming weights)
#   middle layer: ranked by |row j of the linearised map from that layer|
# Stratification potential = Var F_j / E Var(F_j | cell)   (oracle cell masses; bias-corrected).
import numpy as np, sys
rng = np.random.default_rng(0)
d = n = 256
def net(L, seed):
    r = np.random.default_rng(seed)
    return [r.standard_normal((n, n)).astype(np.float32)*np.sqrt(2/n) for _ in range(L)]
def ratio(F, cell, C):
    cnt = np.bincount(cell, minlength=C).astype(float); s1 = np.bincount(cell, F, minlength=C); s2 = np.bincount(cell, F*F, minlength=C)
    ok = cnt > 1
    within = ((s2[ok] - s1[ok]**2/cnt[ok])/(cnt[ok]-1)*cnt[ok]).sum()/cnt[ok].sum()
    return F.var()/within
for L in [int(a) for a in sys.argv[1:]] or [4, 8, 16]:
    Ws = net(L, 200 + L)
    N = 1 << 20
    X = rng.standard_normal((N, d)).astype(np.float32)
    H = X; Z = []
    for W in Ws:
        z = H @ W.T; Z.append(z > 0); H = np.maximum(z, 0)
    F = H.astype(np.float64)
    J = [j for j in range(n) if F[:, j].var() > 1e-8][:8]
    lin = [None]*L                          # linearised map from layer-l hidden units to the output pre-activations
    M = np.eye(n)
    for l in range(L-1, -1, -1):
        lin[l] = M                          # d z_L / d h_l (gates replaced by 1/2), shape (n, n)
        M = 0.5*M @ Ws[l].astype(np.float64) if l > 0 else M
    rL = M @ Ws[0].astype(np.float64)       # unused placeholder
    res = {}
    for name, l in [('layer 1', 0), ('middle', L//2 - 1), ('layer L-1', L-2)]:
        for k in [4, 8, 12]:
            gains = []
            for j in J:
                if name == 'layer 1':
                    rj = (lin[0] @ Ws[0].astype(np.float64))[j] if False else None
                    # linearised row of output j in input space: rows of W_1 weighted by lin map
                    score = np.abs(lin[0][j] @ np.diag(np.ones(n)))          # sensitivity of z_L,j to h_1,i
                elif name == 'layer L-1':
                    score = np.abs(Ws[-1][j].astype(np.float64))
                else:
                    score = np.abs(lin[l][j])
                idx = np.argsort(-score)[:k]
                bits = Z[l][:, idx].astype(np.int64)
                cell = (bits * (1 << np.arange(k))).sum(1)
                gains.append(ratio(F[:, j], cell, 1 << k))
            res[(name, k)] = np.mean(gains)
    print(f"L={L}: " + " | ".join(f"{nm} k={k}: {g:.3f}" for (nm, k), g in res.items()), flush=True)
