# Strand-1 gate test: is the coherent (mixture-shaped) part of the sources' third cumulant concentrated in rank?
# From the live legs of layer L (v29legs_off{net}.pkl, uncompressed n x n legs), form the (2,1) slice of the sources'
# kappa3 tensor, K21_ac = T_aac = sum_s sum_j [w2_j (A_aj^2 P_cj + 2 A_aj P_aj A_cj)/3 + (2 M_aj P_aj P_cj + M_cj P_aj^2)/9],
# project it on the two coherent patterns of the scale mixture / skew (mu_a S_ac and sigma_a^2 mu_c), and compare the
# singular spectrum of K21 with that of the residual K21 - coherent: the ranks needed for 50 / 90 / 99% of the Frobenius
# mass.  If stripping the coherent part lowers the rank materially, carrying it analytically saves transport ranks.
import numpy as np, pickle, sys, time
net = int(sys.argv[1]) if len(sys.argv) > 1 else 0; layer = int(sys.argv[2]) if len(sys.argv) > 2 else 10
P = pickle.load(open(f"v29legs_off{net}.pkl", "rb")); legs = {d["layer"]: d for d in P["legs"]}; dumps = {d["layer"]: d for d in P["dumps"]}
Lg = legs[layer]; d = dumps[layer]; mu, var, C_off = d["mu"].astype(np.float64), d["var"].astype(np.float64), d["C_off"].astype(np.float64); n = len(mu)
S = C_off.copy(); np.fill_diagonal(S, var)
A, Pl, Z, L = Lg["A"].astype(np.float64), Lg["P"].astype(np.float64), Lg["Z"].astype(np.float64), Lg["L"].astype(np.float64)
k = A.shape[0]; w2b = [np.asarray(x, dtype=np.float64) for x in Lg["w2b"]]; s_l = [np.asarray(x, dtype=np.float64) for x in Lg["s"]]; e_l = [np.asarray(x, dtype=np.float64) for x in Lg["e"]]
M = [Pl[s] * s_l[s][None, :] + 3.0 * A[s] * e_l[s][None, :] + Z[s] @ L[s].T for s in range(k)]
def slice21(As, Ps, Ms, w2):
    return ((As * As * w2[None, :]) @ Ps.T + 2.0 * ((As * Ps * w2[None, :]) @ As.T)) / 3.0 + (2.0 * ((Ms * Ps) @ Ps.T) + (Ps * Ps) @ Ms.T) / 9.0
def ranks(Mx):
    sv = np.linalg.svd(Mx, compute_uv=False); c = np.cumsum(sv**2) / np.sum(sv**2)
    return [int(np.searchsorted(c, q) + 1) for q in (0.5, 0.9, 0.99)], sv
f1 = mu[:, None] * S; f2 = var[:, None] * mu[None, :]
def coherent(K):
    X = np.stack([f1.ravel(), f2.ravel()], 1); coef, *_ = np.linalg.lstsq(X, K.ravel(), rcond=None); C = (X @ coef).reshape(n, n)
    return C, coef, 1 - np.sum((K - C)**2) / np.sum(K**2)
t0 = time.time(); K = np.zeros((n, n)); per = []
print(f"net {net} layer {layer}: {k} sources, n = {n}")
print("  source | age | coherent share of its (2,1) slice (R^2 on mu_a S_ac, sigma_a^2 mu_c) | ranks for 50/90/99% mass: slice | residual")
for s in range(k):
    Ks = slice21(A[s], Pl[s], M[s], w2b[s]); K += Ks
    Cs, coef, r2 = coherent(Ks); ra, _ = ranks(Ks); rb, _ = ranks(Ks - Cs)
    print(f"  {s:2d} | {layer - s if s < layer else 0:3d} | {r2:.3f} (g1 {coef[0]:+.4f}, g2 {coef[1]/0.5:+.4f}) | {ra} | {rb}   [{time.time()-t0:.0f}s]", flush=True)
C, coef, r2 = coherent(K); ra, sva = ranks(K); rb, svb = ranks(K - C)
print(f"\n  all sources: coherent share {r2:.3f} (g1 {coef[0]:+.4f}, g2 {coef[1]/0.5:+.4f}; chain D21 fit in note XVIII's record: g ~ 0.02 at depth)")
print(f"  ranks for 50/90/99% of Frobenius mass: full slice {ra} | after stripping the coherent part {rb}")
print("  singular values (full / residual) at ranks 1, 8, 64, 128, 256, 384: " + " ".join(f"{sva[i]:.2e}/{svb[i]:.2e}" for i in (0, 7, 63, 127, 255, 383)))
# the chain's own D21 at this layer for comparison (what it feeds forward), if dumped
if "D21" in Lg:
    D21 = Lg["D21"].astype(np.float64); Cd, coefd, r2d = coherent(D21); rd, _ = ranks(D21); rdr, _ = ranks(D21 - Cd)
    print(f"  chain's carried D21: coherent share {r2d:.3f}, ranks {rd} -> residual {rdr}; corr(D21, K21 from legs) {np.corrcoef(D21.ravel(), K.ravel())[0,1]:+.3f}")
