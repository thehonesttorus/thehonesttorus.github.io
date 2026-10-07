# Exact check (note XXXVI): slices of the linear transport W# of a pure (3,1)-class tensor S31(B), B_ab = T_aaab (a != b),
# polynomial 4 sum_(a != b) B_ab v_a^3 v_b. With G = B W^T (G_aj = sum_b B_ab W_jb), H = W o W, W3 = W o H:
#   diag_i   = 4 sum_a W_ia^3 G_ai
#   (2,2)_ij = 2 [(H X^T) + (H X^T)^T]_ij,   X = W o G^T   (X_ja = W_ja G_aj)
#   (3,1)_ij = kappa(i,i,i,j) = (W3 G)_ij + 3 ((H o G^T) W^T)_ij
import itertools, numpy as np
rng = np.random.default_rng(11); m, n = 6, 5
W = rng.standard_normal((n, m)); B = rng.standard_normal((m, m)); np.fill_diagonal(B, 0.0)
T = np.zeros((m,) * 4)
for a, b in itertools.permutations(range(m), 2):
    for idx in set(itertools.permutations((a, a, a, b))):
        T[idx] = B[a, b]
Tz = np.einsum("ia,jb,kc,ld,abcd->ijkl", W, W, W, W, T)
H = W * W; W3 = W * H; G = B @ W.T; X = W * G.T
diag = 4 * np.einsum("ia,ai->i", W3, G)
K22 = 2 * (H @ X.T + (H @ X.T).T)
K31 = W3 @ G + 3 * (H * G.T) @ W.T
off = ~np.eye(n, dtype=bool)
print("diag", np.abs(diag - np.array([Tz[i, i, i, i] for i in range(n)])).max(),
      "| (2,2)", np.abs((K22 - np.array([[Tz[i, i, j, j] for j in range(n)] for i in range(n)]))[off]).max(),
      "| (3,1)", np.abs((K31 - np.array([[Tz[i, i, i, j] for j in range(n)] for i in range(n)]))[off]).max())
