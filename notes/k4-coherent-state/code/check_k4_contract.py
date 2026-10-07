# Exact checks for note XXXII (no network, no Monte Carlo): the fourth-cumulant contract of est_v29 (V33_K4Q = 3, the
# regenerated core with adaptive lambda), written as tensors.  python check_k4_contract.py
#   (1) the code's diagonal  g4row = METRIC_C * (dG + k4corr)  equals  Q(K, d) - R_W(E) + 2 lam s_off^2,
#       R_W(E) = 3 diag(H E H^T) - 2 H P(E 1),  E = K - Khat0  (Khat = the Ritz approximation, zero diagonal),
#       and THEORY_3's split R_W(E) = D_metric + D_known + D_unresolved sums back to it;
#   (2) every slice of the exactly transported retained tensor W#T, T = diag class (d) + S22(K):
#       diag = Q = (W^o4) d + 3 diag(H K H^T);  (2,2)_ij = (H K H^T)_ij + (H D H^T)_ij + 2 u_ij^T K u_ij,  u_ij = W_i o W_j;
#       (3,1)_ij = kappa(i,i,i,j) = [((H A) o W) W^T]_ij,  A = 3K + diag(d)   -- against an explicit 4-tensor;
#   (3) wk4m = (dG_i + dG_j)/3 is the (2,2) slice of the core J(2I, N) with diag N = dG, whose (3,1) slice is N_ij;
#       the code's wk431 = lam C_off replaces the pair core's part of N_ij.
import itertools, numpy as np
rng = np.random.default_rng(7)
m, n = 6, 5                                   # input (y) width m, output (z) width n; the code has m = n
W = rng.standard_normal((n, m)) * np.sqrt(2 / m)
K = rng.standard_normal((m, m)); K = K + K.T; np.fill_diagonal(K, 0.0)         # y-level (2,2) slice, zero diagonal
d = rng.standard_normal(m)                                                      # y-level kappa_4 diagonal
H = W * W; W4 = H * H; one = np.ones(m)
cA, cI = 6.0 / (m + 4.0), -3.0 / ((m + 2.0) * (m + 4.0))                          # st["cA"], st["cI"] at width m
P = lambda v: cA * v + cI * np.sum(v) * one
# --- (1) the code path, operation by operation (est_v29 regen block, K4Q = 3, BETA != 0) ---
V = np.linalg.qr(rng.standard_normal((m, 3)))[0]; lb = rng.standard_normal(3)   # any Ritz pairs (code: 12 of them)
k22q, k4q = K, d
g_prev = (k4q + k22q @ one) * cA + (np.sum(k4q) + np.sum(k22q @ one)) * cI
k4corr = (W4 @ k4q) * 0.5 - H @ (k4q * cA + np.sum(k4q) * cI)
WV = H @ V; Kd = (V * V) @ lb; rs = V @ (lb * V.sum(0)) - Kd
k4corr = k4corr + ((WV * WV) @ lb) * 1.5 - (W4 @ Kd) * 1.5 - H @ (rs * cA + np.sum(rs) * cI)
lam = 0.37; var_prev = rng.random(m) + 0.5; Cy = np.diag(var_prev) + 0.1 * (lambda X: X + X.T)(rng.standard_normal((m, m)))
np.fill_diagonal(Cy, var_prev); var = np.diag(W @ Cy @ W.T)
t_g = H @ g_prev; t_v = var - H @ var_prev; dG = t_g + t_v * lam
g4row = dG * 2.0 + k4corr * 2.0
Q = W4 @ d + 3 * np.diag(H @ K @ H.T)
Khat0 = V @ np.diag(lb) @ V.T; np.fill_diagonal(Khat0, 0.0); E = K - Khat0
RW = lambda E_: 3 * np.diag(H @ E_ @ H.T) - 2 * H @ P(E_ @ one)
print("(1) g4row - [Q - R_W(E) + 2 lam s_off^2] =", np.abs(g4row - (Q - RW(E) + 2 * lam * t_v)).max())
s = E @ one; S = s.sum(); g = (s - S / (2 * (m - 1)) * one) / (m - 2)
Ag = np.outer(g, one) + np.outer(one, g) - 2 * np.diag(g); E0 = E - Ag; r = H @ one; p = P(s)
Dm = (r - 2) * (H @ p); Dk = 3 * np.diag(H @ Ag @ H.T) - r * (H @ p); Du = 3 * np.diag(H @ E0 @ H.T)
print("    row-null E0 1 =", np.abs(E0 @ one).max(), "| D_metric + D_known + D_unresolved - R_W(E) =", np.abs(Dm + Dk + Du - RW(E)).max())
# --- (2) explicit tensor T = sum_a d_a e_a^4 + S22(K), transported by W on all four legs ---
T = np.zeros((m,) * 4)
for a in range(m):
    T[a, a, a, a] = d[a]
for a, b in itertools.permutations(range(m), 2):
    for idx in set(itertools.permutations((a, a, b, b))):
        T[idx] = K[a, b]
Tz = np.einsum("ia,jb,kc,ld,abcd->ijkl", W, W, W, W, T)
dz = np.array([Tz[i, i, i, i] for i in range(n)])
K22z = np.array([[Tz[i, i, j, j] for j in range(n)] for i in range(n)])
K31z = np.array([[Tz[i, i, i, j] for j in range(n)] for i in range(n)])
off = ~np.eye(n, dtype=bool); A = 3 * K + np.diag(d)
K22f = H @ K @ H.T + H @ np.diag(d) @ H.T + 2 * np.einsum("ia,ja,ab,ib,jb->ij", W, W, K, W, W)
K31f = ((H @ A) * W) @ W.T
print("(2) diag vs Q:", np.abs(dz - Q).max(), "| (2,2):", np.abs((K22z - K22f)[off]).max(), "| (3,1):", np.abs((K31z - K31f)[off]).max())
# --- (3) wk4m as the (2,2) slice of J(2I, N): K_ij = (M_ii N_jj + M_jj N_ii + 4 M_ij N_ij)/6 with M = 2I ---
wk4m = (dG[:, None] + dG[None, :]) / 3.0
print("(3) wk4m - (2,2) slice of J(2I, N) with diag N = dG:", np.abs((wk4m - (2 * dG[None, :] + 2 * dG[:, None]) / 6)[off]).max())
