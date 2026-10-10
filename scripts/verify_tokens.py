"""Identity checks for the token (Wiener-chaos / Fock) picture, notes/stage10/s9_tokens.tex. No estimator runs.
  python scripts/verify_tokens.py"""
import numpy as np
from math import sqrt, pi, acos, cos
from scipy.stats import norm
rng = np.random.default_rng(11)
J = lambda th: np.sin(th) + (pi - th) * np.cos(th)

print("A. two-layer network: factorised kink current vs Gaussian closure vs truth; the contact term")
n, N = 256, 1 << 20
W1 = rng.standard_normal((n, n)) * sqrt(2 / n); W2 = rng.standard_normal((n, n)) * sqrt(2 / n)
nw = np.linalg.norm(W1, axis=1); G = W1 @ W1.T; cosab = np.clip(G / np.outer(nw, nw), -1, 1); th = np.arccos(cosab)
mu1 = nw / sqrt(2 * pi)                                    # E h_1a (exact)
E11 = np.outer(nw, nw) * J(th) / (2 * pi)                  # E h_a h_b (exact arc-cosine kernel)
C1 = E11 - np.outer(mu1, mu1)
m2 = W2 @ mu1; s2 = np.sqrt(np.einsum('ca,ab,cb->c', W2, C1, W2)); al = m2 / s2
u_gc = m2 * norm.cdf(al) + s2 * norm.pdf(al)               # Gaussian closure for h_2
K0 = (pi - th) / (2 * pi)                                  # P(g_a g_b = 1), exact at layer 1
grad2 = np.einsum('ca,ab,cb->c', W2, K0 * G, W2)           # E |grad z_2c|^2 (exact)
f0 = norm.pdf(al) / s2
u_mf = m2 * norm.cdf(al) + f0 * grad2                      # independent walkers (factorised kink current)
contact = f0 * (W2 ** 2 @ mu1 ** 2)                        # f_c * sum_a W_ca^2 mu_a^2
acc = np.zeros(n)
for _ in range(8):
    x = rng.standard_normal((N // 8, n)); acc += np.maximum(np.maximum(x @ W1.T, 0) @ W2.T, 0).sum(0)
u = acc / N
print("   rms error: Gaussian closure %.2e   independent walkers %.2e   (rms truth %.3f, MC noise ~%.0e)" % (np.sqrt(np.mean((u_gc - u) ** 2)), np.sqrt(np.mean((u_mf - u) ** 2)), np.sqrt(np.mean(u ** 2)), 1.4 / sqrt(N)))
print("   E|grad z|^2 - sigma^2: mean %.3f ; contact sum sum_a W_ca^2 mu_a^2: mean %.3f ; corr %.3f" % (np.mean(grad2 - s2 ** 2), np.mean(W2 ** 2 @ mu1 ** 2), np.corrcoef(grad2 - s2 ** 2, W2 ** 2 @ mu1 ** 2)[0, 1]))
print("   (u_mf - u_gc) vs contact term f_c sum W^2 mu^2: slope %.3f corr %.3f ; after removing it: rms vs GC %.2e" % (np.polyfit(contact, u_mf - u_gc, 1)[0], np.corrcoef(contact, u_mf - u_gc)[0, 1], np.sqrt(np.mean((u_mf - contact - u_gc) ** 2))))

print("B. layer-2 cumulants are walk sums on the hidden graph G = W1 W1^T (chaos <= 2 truncation), n = 24")
n, N = 24, 1 << 23
W1 = rng.standard_normal((n, n)) * sqrt(2 / n); W2 = rng.standard_normal((n, n)) * sqrt(2 / n); G = W1 @ W1.T
nw = np.linalg.norm(W1, axis=1); c = 0; om = W2[c]; d = om / (2 * sqrt(2 * pi) * nw); Q = W1.T @ (d[:, None] * W1); f1 = 0.5 * W1.T @ om
k3_walk = 6 * f1 @ Q @ f1 + 8 * np.trace(Q @ Q @ Q)
k3_hub = 6 * 0.25 * om @ G @ (d[:, None] * G) @ om            # same open walk written on hidden neurons
mom = np.zeros(3)
for _ in range(16):
    x = rng.standard_normal((N // 16, n)); z = np.maximum(x @ W1.T, 0) @ om
    mom += [z.sum(), (z ** 2).sum(), (z ** 3).sum()]
mom /= N; k3_mc = mom[2] - 3 * mom[1] * mom[0] + 2 * mom[0] ** 3
print("   kappa_3: walk formula %.5f (open part on hidden graph %.5f = %.5f), Monte Carlo %.5f" % (k3_walk, 6 * f1 @ Q @ f1, k3_hub, k3_mc))
print("   kappa_2: walk formula |f1|^2 + 2 tr Q^2 = %.5f (+ chaos>=4: %.5f), Monte Carlo %.5f" % (f1 @ f1 + 2 * np.trace(Q @ Q), (om ** 2 @ nw ** 2) * (0.5 - 1 / (2 * pi) - 0.25 - 1 / (4 * pi)) * 0 + (om @ ((np.outer(nw, nw) * J(np.arccos(np.clip(G / np.outer(nw, nw), -1, 1))) / (2 * pi) - np.outer(nw, nw) / (2 * pi)) @ om)) - f1 @ f1 - 2 * np.trace(Q @ Q), mom[1] - mom[0] ** 2))

print("C. mean token number of pre-activations, N = E|grad z|^2 / Var z: infinite width 1/(1 - cos theta) vs a width-256 He net")
th, pred = pi / 2, [1.0]
for l in range(1, 16): th = acos(min(1.0, J(th) / pi)); pred.append(1 / (1 - cos(th)))
n, L, S, R = 256, 16, 1 << 14, 96
Ws = [rng.standard_normal((n, n)) * sqrt(2 / n) for _ in range(L)]
x = rng.standard_normal((S, n)); h = x; var = []
for W in Ws:
    z = h @ W.T; var.append(z.var(0).mean()); h = np.maximum(z, 0)
g2 = np.zeros(L)
for s in range(R):
    Jm = Ws[0].copy(); hh = Ws[0] @ x[s]; g2[0] += np.mean(np.sum(Jm ** 2, 1))
    for l in range(1, L):
        Jm = Ws[l] @ ((hh > 0)[:, None] * Jm); hh = Ws[l] @ np.maximum(hh, 0); g2[l] += np.mean(np.sum(Jm ** 2, 1))
g2 /= R
for l in (0, 1, 2, 4, 8, 12, 15):
    print("   layer %2d: infinite-width %.2f   width-256 net %.2f" % (l + 1, pred[l], g2[l] / var[l]))
