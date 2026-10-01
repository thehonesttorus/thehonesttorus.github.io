"""E0b: exactness of I2 with a finite-variance estimator (replaces the line method of e0, whose 1/|v.grad z|
weights have infinite variance).  delta(z) -> Gaussian kernel of width h*s_km, two widths, Richardson
(bias O(h), since H jumps where other walls cross; linear extrapolation).  RHS_lj = sum_{k<=l, m} E[K_h(z_km) |grad z_km|^2 S_(l,j)<-(k,m)];  LHS = E a_lj.
Paired on the same samples; standard errors from batch means.  Usage: seed n L N"""
import numpy as np, sys
seed, n, L, N = [int(float(a)) for a in (sys.argv[1:] + ['0', '8', '3', '4e6'][len(sys.argv) - 1:])][:4]
rng = np.random.default_rng(seed)
W = [rng.standard_normal((n, n)) * np.sqrt(2 / n) for _ in range(L)]
X = rng.standard_normal((200000, n)); s = []
a = X
for l in range(L):
    z = a @ W[l]; s.append(z.std(0)); a = np.maximum(z, 0)
hs = (0.01, 0.02)
B = 20000; nb = N // B
lhs_b = np.zeros((nb, L, n)); rhs_b = np.zeros((2, nb, L, n))
for b in range(nb):
    X = rng.standard_normal((B, n)); a = X; J = np.broadcast_to(np.eye(n), (B, n, n)).copy()
    zs, Gs, gs = [], [], []
    for l in range(L):
        z = a @ W[l]; G = J @ W[l]; g = (z > 0).astype(float)
        zs.append(z); Gs.append(G); gs.append(g); a = z * g; J = G * g[:, None, :]
        lhs_b[b, l] = a.mean(0)
    for k in range(L):
        g2 = (Gs[k] ** 2).sum(1)                                   # (B, n): |grad z_km|^2
        for hi, h in enumerate(hs):
            bw = h * s[k]
            base = np.exp(-0.5 * (zs[k] / bw) ** 2) / (np.sqrt(2 * np.pi) * bw) * g2   # (B, m)
            rhs_b[hi, b, k] += base.mean(0)
            # downstream: sens (B, m, j) = d a_l / d a_km
            sens = np.broadcast_to(np.eye(n), (B, n, n))
            for l in range(k + 1, L):
                sens = (sens @ W[l]) * gs[l][:, None, :]
                rhs_b[hi, b, l] += np.einsum('bm,bmj->j', base, sens) / B
lhs = lhs_b.mean(0); r1, r2 = rhs_b.mean(1)
rich = 2 * r1 - r2   # bias is O(h) (H jumps at crossing walls): linear extrapolation
d = lhs_b - (2 * rhs_b[0] - rhs_b[1])
se = d.std(0, ddof=1) / np.sqrt(nb)
print(f"n={n} L={L} N={nb*B}; linear extrapolation over kernel widths {hs} x s")
print("layer-1 closed form check: max|LHS-|w|/sqrt(2pi)| = %.1e" % np.abs(lhs[0] - np.linalg.norm(W[0], axis=0) / np.sqrt(2 * np.pi)).max())
for l in range(L):
    print("layer %d: max|LHS-RHS| %.1e (rel %.1e)  max |z| %.2f   [unextrapolated h=%.2f: %.1e]" % (
        l + 1, np.abs(lhs[l] - rich[l]).max(), np.abs(lhs[l] - rich[l]).max() / np.abs(lhs[l]).max(),
        (np.abs(lhs[l] - rich[l]) / se[l]).max(), hs[0], np.abs(lhs[l] - r1[l]).max()))
