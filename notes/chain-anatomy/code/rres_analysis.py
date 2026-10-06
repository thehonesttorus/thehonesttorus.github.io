# What is in the slice residual Rres = K21 - rep21 - S_sep that the chain truncates to rank 16?
# Hypothesis: the mismatch between the legs' Phi^2 Phi transport of the carried (2,1) slice and the exact pair-level
# coefficient, i.e. a diagonal rescaling d(x) D21 d(y) of the carried slice (plus Gaussian Mehler-2 terms).
import numpy as np, pickle, sys
from scipy.special import ndtr
net = int(sys.argv[1]) if len(sys.argv) > 1 else 0
D = pickle.load(open(f"v29dump_off{net}.pkl", "rb"))
phi = lambda x: np.exp(-0.5 * x * x) / np.sqrt(2 * np.pi)
def energy(M, ranks):
    s = np.linalg.svd(M, compute_uv=False); e = np.cumsum(s**2) / np.sum(s**2)
    return s, [e[r - 1] for r in ranks]
def bilinear_fit(target, M, feats, names):
    """least squares target ~ sum_pq beta_pq d(f_p) M d(f_q) (off-diagonal entries)"""
    n = M.shape[0]; mask = ~np.eye(n, dtype=bool)
    cols = []; labels = []
    for p, fp in enumerate(feats):
        for q, fq in enumerate(feats):
            cols.append(((fp[:, None] * M) * fq[None, :])[mask]); labels.append(f"{names[p]}*{names[q]}")
    A = np.stack(cols, 1); y = target[mask]
    beta, *_ = np.linalg.lstsq(A, y, rcond=None); pred = A @ beta
    return 1 - np.sum((y - pred)**2) / np.sum(y**2), beta, labels, pred, mask
for d in D:
    l = d["layer"]; n = len(d["mu"])
    mu, var = d["mu"], d["var"]; sig = np.sqrt(var); a = mu / sig; Phi = ndtr(a); ph = phi(a); m = sig * (a * Phi + ph)
    Rres, K21, D21, D21w, S_sep, C_off = d["Rres"], d["K21"], d["D21"], d["D21_w"], d["S_sep"], d["C_off"]
    s_r, e_r = energy(Rres, [16, 64, 256]); s_k, e_k = energy(K21, [16, 64, 256]); s_d, e_d = energy(D21, [16, 64, 256])
    print(f"\nlayer {l}: |Rres|_F {np.linalg.norm(Rres):.3e}  |K21| {np.linalg.norm(K21):.3e}  |rep21| {np.linalg.norm(d['rep21']):.3e}  |S_sep| {np.linalg.norm(S_sep):.3e}  |D21| {np.linalg.norm(D21):.3e}")
    print(f"   energy fraction in top 16/64/256 singular directions: Rres {e_r[0]:.3f}/{e_r[1]:.3f}/{e_r[2]:.3f}   K21 {e_k[0]:.3f}/{e_k[1]:.3f}/{e_k[2]:.3f}   D21 {e_d[0]:.3f}/{e_d[1]:.3f}/{e_d[2]:.3f}")
    # candidate 1: diagonal rescalings of the carried slice D21 (and its transpose)
    feats = [np.ones(n), Phi, ph, Phi * Phi, ph * m / sig, m / sig, a * ph]; names = ["1", "Phi", "phi", "Phi2", "phim/s", "m/s", "aphi"]
    r2a, beta, labels, pred_a, mask = bilinear_fit(Rres, D21, feats, names)
    r2b, _, _, pred_b, _ = bilinear_fit(Rres, D21.T, feats, names)
    # candidate 2: Mehler-2 Gaussian term d(x) (C_off*C_off) d(y)
    r2c, _, _, pred_c, _ = bilinear_fit(Rres, C_off * C_off, feats, names)
    # candidate 3: the gain-form mean coupling d(x) (mu (x) S-ish: mu_i C_off_ij) d(y) and sigma^2_i mu_j
    r2d, _, _, pred_d, _ = bilinear_fit(Rres, mu[:, None] * C_off, feats, names)
    # joint fit of all four families
    cols = []
    for M in (D21, D21.T, C_off * C_off, mu[:, None] * C_off, var[:, None] * mu[None, :] * (1 - np.eye(n))):
        for fp in feats:
            for fq in feats: cols.append(((fp[:, None] * M) * fq[None, :])[mask])
    A = np.stack(cols, 1); y = Rres[mask]; b, *_ = np.linalg.lstsq(A, y, rcond=None); pj = A @ b; r2j = 1 - np.sum((y - pj)**2) / np.sum(y**2)
    print(f"   R^2 of Rres against: d(.) D21 d(.) {r2a:.3f} | d(.) D21^T d(.) {r2b:.3f} | d(.) (C_off*C_off) d(.) {r2c:.3f} | d(.) mu_i C_ij d(.) {r2d:.3f} | all jointly {r2j:.3f}")
    # the explicit prediction: exact pair coefficient minus the legs' Phi^2: c_i = Phi_i - phi_i m_i/sig_i - Phi_i^2, times Phi_j
    c = Phi - ph * m / sig - Phi * Phi
    pred_exact = (c[:, None] * D21) * Phi[None, :]; np.fill_diagonal(pred_exact, 0)
    coef = np.sum(pred_exact * Rres) / np.sum(pred_exact**2); r2e = 1 - np.sum((Rres - coef * pred_exact)**2) / np.sum(Rres**2)
    print(f"   explicit mismatch d(Phi - phi m/s - Phi^2) D21 d(Phi): best scale {coef:+.3f}, R^2 {r2e:.3f}; same with D21^T: ", end="")
    pe2 = (c[:, None] * D21.T) * Phi[None, :]; np.fill_diagonal(pe2, 0); coef2 = np.sum(pe2 * Rres) / np.sum(pe2**2); print(f"scale {coef2:+.3f}, R^2 {1 - np.sum((Rres - coef2*pe2)**2)/np.sum(Rres**2):.3f}")
    # what the rank-16 truncation keeps vs drops, and whether the dropped part is the explained (physical) part
    U, s, Vt = np.linalg.svd(Rres); R16 = (U[:, :16] * s[:16]) @ Vt[:16]; drop = Rres - R16
    full_pred = np.zeros((n, n)); full_pred[mask] = pj
    print(f"   joint fit explains: of the top-16 part {1 - np.sum((R16 - full_pred)[mask]**2)/np.sum(R16[mask]**2):.3f} (by overlap), of the dropped part {np.sum(drop[mask]*full_pred[mask])/np.sum(drop[mask]**2):.3f} (regression slope); |drop|/|Rres| {np.linalg.norm(drop)/np.linalg.norm(Rres):.3f}")
    # alignment of Rres's top left/right singular vectors with the mean direction and the top eigenvector of C_off
    ev, V = np.linalg.eigh(C_off + np.diag(var)); u1 = V[:, -1]; mhat = mu / np.linalg.norm(mu)
    print(f"   top singular vectors of Rres: |<u1,mean>| left {abs(U[:,0]@mhat):.3f} right {abs(Vt[0]@mhat):.3f}; |<.,topeig C>| left {abs(U[:,0]@u1):.3f} right {abs(Vt[0]@u1):.3f}; singular values s1..s4/s17 {s[:4]/s[16]}")
