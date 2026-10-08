# Two measurements behind the deck (Moreau) symmetry note: per layer of one network's free-running chain dump,
#   (1) the alpha = mu / sigma distribution: the fraction of units saturated off (alpha <= -t), saturated on
#       (alpha >= t) and marginal, for t = 2, 2.5, 3; relu(z) = z + relu(-z) maps an on-saturated unit to an
#       off-saturated one with the same hinge coefficients phi(alpha)/sigma;
#   (2) how much of the (2,1) third-cumulant slice D21[a, b] = kappa(z_a, z_a, z_b) the pair program can read: its
#       rows enter with the hinge coefficient c_a(1,2) = phi(alpha_a)/sigma_a, so the share of the read-weighted
#       energy sum_ab (c_a D21_ab)^2 carried by rows with |alpha_a| >= t is the error of dropping those rows;
#   (3) the singular spectrum of D21 and of its read-weighted form diag(c) D21 diag(Phi): the effective rank a
#       randomized sketch of the readout would need (energy captured at ranks 32 ... 512).
#   python deckspec.py DUMP.npz NET
import sys, numpy as np
from math import erf, sqrt, pi

d = np.load(sys.argv[1]); net = int(sys.argv[2])
Wc = np.load(f"../official/W_off{net}.npy").astype(np.float64)   # column convention: z_l = Wc[l] @ y_(l-1)
# the dump's "mu" at layers < 15 is the post-activation mean m_l (written after the nonlinearity); the pre-activation
# mean is mu_l = Wc[l] @ m_(l-1); the last layer's dump is taken before the nonlinearity (checked below)
mpost = {l: d[f"mu_{l}"].astype(np.float64) for l in range(15) if f"mu_{l}" in d.files}
def premean(l):
    if l == 0:
        return np.zeros(Wc.shape[1])
    return Wc[l] @ mpost[l - 1]
if "mu_15" in d.files and 14 in mpost:
    print(f"check: last-layer dump mu vs Wc[15] @ m_14: rel diff {np.linalg.norm(d['mu_15'] - Wc[15] @ mpost[14]) / np.linalg.norm(d['mu_15']):.2e}")
phi = lambda a: np.exp(-0.5 * a * a) / sqrt(2 * pi)
Phi = np.vectorize(lambda a: 0.5 * (1 + erf(a / sqrt(2))))
print(f"net {net}: dump keys sample {sorted(d.files)[:12]}")
for l in range(16):
    if f"mu_{l}" not in d.files or f"var_{l}" not in d.files:
        continue
    mu = premean(l) if l < 15 else d[f"mu_{l}"].astype(np.float64); var = d[f"var_{l}"].astype(np.float64)
    s = np.sqrt(np.maximum(var, 1e-12)); al = mu / s
    fr = " ".join(f"t={t}: off {np.mean(al <= -t):.3f} on {np.mean(al >= t):.3f}" for t in (2.0, 2.5, 3.0))
    line = f"layer {l:2d}  alpha sd {al.std():.2f}  {fr}"
    if f"D21_{l}" in d.files:
        D = d[f"D21_{l}"].astype(np.float64)
        c = phi(al) / s; P = Phi(al)
        R = c[:, None] * D * P[None, :]                     # what the pair program reads (row-weighted by the hinge)
        e = (R * R).sum(axis=1); tot = e.sum()
        drop = " ".join(f"|a|>={t}: {e[np.abs(al) >= t].sum() / tot:.2e}" for t in (2.0, 2.5, 3.0))
        sv = np.linalg.svd(D, compute_uv=False) ** 2; svr = np.linalg.svd(R, compute_uv=False) ** 2
        cap = lambda v: " ".join(f"{k}:{v[:k].sum() / v.sum():.4f}" for k in (32, 64, 128, 256, 512))
        line += f"\n     read-energy share of dropped rows {drop}\n     D21 energy at rank {cap(sv)}\n     read-weighted at rank {cap(svr)}"
    print(line, flush=True)

# (4) is the read-weighted slice's dominant subspace one the chain already has?  Right subspace: the top-r eigenvectors
# of the gated covariance Phi C Phi (C = C_off + diag var, the next transport's input); left: of c C c. Energy of
# R captured by R U U^T (right), U' U'^T R (left) and U' U'^T R U U^T (both), against the optimal rank-r SVD value.
print("alignment of the read-weighted slice with the gated covariance eigenvectors (captured energy, r = 32/64/128/256)")
for l in range(1, 15):
    if f"D21_{l}" not in d.files or f"C_off_{l}" not in d.files:
        continue
    mu = premean(l); var = d[f"var_{l}"].astype(np.float64); s = np.sqrt(np.maximum(var, 1e-12)); al = mu / s
    c = phi(al) / s; P = Phi(al)
    R = c[:, None] * d[f"D21_{l}"].astype(np.float64) * P[None, :]
    C = d[f"C_off_{l}"].astype(np.float64); C = C - np.diag(np.diag(C)) + np.diag(var)
    er, Ur = np.linalg.eigh(P[:, None] * C * P[None, :]); Ur = Ur[:, ::-1]
    el, Ul = np.linalg.eigh(c[:, None] * C * c[None, :]); Ul = Ul[:, ::-1]
    sv = np.linalg.svd(R, compute_uv=False) ** 2; tot = sv.sum(); out = []
    for r in (32, 64, 128, 256):
        U = Ur[:, :r]; V = Ul[:, :r]
        right = np.sum((R @ U) ** 2) / tot; left = np.sum((V.T @ R) ** 2) / tot; both = np.sum((V.T @ R @ U) ** 2) / tot
        out.append(f"r{r}: svd {sv[:r].sum() / tot:.4f} right {right:.4f} left {left:.4f} both {both:.4f}")
    print(f"layer {l:2d}  " + "  ".join(out), flush=True)
