"""Escape hatch 1: can the old third-order tensor T_t = sum_r w_r y_r (x) y_r (x) z_r (age >= 3, all sources aggregated)
be carried by k aggregated signed n x n states A_j = T_t[., ., u_j]?  Spectrum of its third-leg Gram
G = Z^T diag(w) (Y Y^T)^{o2} diag(w) Z: fraction of ||T||_F^2 in the top-k third-leg modes, and the D21 relative error
of the rank-k (third-leg) truncation. usage: python t_hatch.py SET MLP targets"""
import sys, json
import numpy as np
sys.path.insert(0, "../../fresh-slate/bench")
import bench, fc_hook as F
name, mlp = sys.argv[1], int(sys.argv[2]); targets = [int(x) for x in sys.argv[3].split(",")]
S = bench.load_set(name); W = bench.weights(S, mlp)
def hook(t, sources, D21):
    if t not in targets: return False
    old = [s for s in sources if t - s["s"] >= 3]; n = D21.shape[0]
    Y = np.concatenate([s["SP"].astype(np.float32) @ s["Z"] for s in old]).astype(np.float64)
    Z = np.concatenate([s["Z"] for s in old]).astype(np.float64); w = np.concatenate([s["w2"] for s in old]).astype(np.float64)
    K = (Y @ Y.T) ** 2 * np.outer(w, w)
    G = Z.T @ K @ Z; ev, U = np.linalg.eigh(G); ev, U = ev[::-1], U[:, ::-1]
    frac = {k: float(ev[:k].sum() / ev.sum()) for k in (1, 8, 64, n // 8, n // 4)}
    # D21 first term only (Y^2 Z) for the truncation test: D = (Y o Y w)^T Z ; truncated: Z -> Z U_k U_k^T
    Dt = ((Y * Y) * w[:, None]).T @ Z
    err = {}
    for k in (1, 8, 64, n // 8, n // 4):
        Uk = U[:, :k]; err[k] = float(np.linalg.norm(Dt - Dt @ Uk @ Uk.T) / np.linalg.norm(Dt))
    r = dict(set=name, mlp=mlp, t=t, atoms=len(w), energy_frac_topk=frac, d21_term1_relerr_topk=err)
    print(json.dumps(r), flush=True); return t >= max(targets)
F.HOOK = hook; F.run(W, k4mf=True)
