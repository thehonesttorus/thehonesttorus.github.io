"""Numerical checks for the note `notes/finite-resolution-question-algebra.md`.

These are identity checks on small synthetic bias-free ReLU MLPs with He-Gaussian
weights and standard-Gaussian inputs (the WhestBench setting).  They are not
benchmark results and do not use the official dataset or harness.

Run:  python3 finite_resolution_checks.py
"""
import numpy as np
from scipy import integrate
from scipy.stats import norm
from numpy.polynomial.hermite_e import hermeval

rng = np.random.default_rng(0)
SQ2PI = np.sqrt(2 * np.pi)


def he(k, x):
    c = np.zeros(k + 1); c[k] = 1.0
    return hermeval(x, c)


# ---------------------------------------------------------------------------
# (a) Hermite coefficients of ReLU with respect to N(0,1):
#     c_k = E[z_+ He_k(z)];  claim: c_0 = 1/sqrt(2pi), c_1 = 1/2,
#     c_k = He_{k-2}(0)/sqrt(2pi) for k >= 2 (all even k nonzero, odd k >= 3 zero).
#     Consequence: no finite Wiener-chaos truncation contains a ReLU'd question.
# ---------------------------------------------------------------------------
print("(a) Hermite coefficients of ReLU  c_k = E[z_+ He_k(z)]")
for k in range(0, 11):
    num, _ = integrate.quad(lambda z: z * he(k, z) * np.exp(-z * z / 2) / SQ2PI, 0, np.inf, limit=200)
    pred = 1 / SQ2PI if k == 0 else (0.5 if k == 1 else he(k - 2, 0.0) / SQ2PI)
    print(f"   k={k:2d}  numeric={num:+.6f}  closed form={pred:+.6f}")

# ---------------------------------------------------------------------------
# (b) The resolution-1 effect of a half-space bit.  In the orthonormal basis
#     (1, x_1, ..., x_n) of Wiener chaos of degree <= 1 in L^2(gamma), the
#     compression e_w = P_1 M_{1{w.x>0}} P_1 has the closed form
#        e_w = [[1/2, w^T/(|w| sqrt(2pi))], [w/(|w| sqrt(2pi)), I/2]],
#     spectrum {1/2 - 1/sqrt(2pi), 1/2 (n-1 times), 1/2 + 1/sqrt(2pi)},
#     unsharpness ||e - e^2|| = ||[P_1, M_bit]||^2 = 1/4, and the mean
#     activation is the vacuum-to-weight matrix element <xi_w, e_w Omega>.
# ---------------------------------------------------------------------------
n = 6
w = rng.normal(size=n); wn = w / np.linalg.norm(w)
N = 3_000_000
X = rng.normal(size=(N, n))
bit = (X @ wn > 0).astype(float)
B = np.concatenate([np.ones((N, 1)), X], axis=1)
E = (B.T * bit) @ B / N
pred = np.zeros((n + 1, n + 1)); pred[0, 0] = 0.5
pred[0, 1:] = wn / SQ2PI; pred[1:, 0] = wn / SQ2PI; pred[1:, 1:] = 0.5 * np.eye(n)
print("\n(b) resolution-1 effect of 1{w.x>0}")
print("   max |Monte Carlo - closed form| =", f"{np.abs(E - pred).max():.2e}")
print("   eigenvalues:", np.round(np.sort(np.linalg.eigvalsh(pred)), 4),
      " predicted extremes", round(0.5 - 1 / SQ2PI, 4), round(0.5 + 1 / SQ2PI, 4))
print("   ||e - e^2|| =", f"{np.linalg.norm(pred - pred @ pred, 2):.6f}", "(claim 0.25)")
xi = np.concatenate([[0.0], w])
print("   <xi_w, e_w Omega> =", f"{xi @ pred[:, 0]:.6f}",
      "  E[(w.x)_+] (MC) =", f"{np.mean(np.maximum(X @ w, 0)):.6f}",
      "  |w|/sqrt(2pi) =", f"{np.linalg.norm(w) / SQ2PI:.6f}")

# ---------------------------------------------------------------------------
# (c) Compressed code effects of a small bias-free ReLU net: they are PSD,
#     they sum to P_1 (a POVM on C (+) input space), they do not commute,
#     and the final-layer means are  sum_c A_c E[x 1_c]  (cell Jacobian times
#     the coherence row of the cell's effect).
# ---------------------------------------------------------------------------
n, L = 4, 3
Ws = [rng.normal(size=(n, n)) * np.sqrt(2 / n) for _ in range(L)]
N = 2_000_000
X = rng.normal(size=(N, n)); h = X.copy(); codes = []
for W in Ws:
    z = h @ W.T; codes.append((z > 0).astype(np.int8)); h = np.maximum(z, 0)
code = np.concatenate(codes, axis=1)
keys = code @ (1 << np.arange(code.shape[1]))
uniq, inv = np.unique(keys, return_inverse=True)
B = np.concatenate([np.ones((N, 1)), X], axis=1)
effects = [(B[inv == ci].T @ B[inv == ci]) / N for ci in range(len(uniq))]
S = sum(effects)
print(f"\n(c) net n={n}, L={L}: {len(uniq)} realised joint codes")
print("   max |sum_c e_c - P_1| =", f"{np.abs(S - np.eye(n + 1)).max():.2e}")
print("   min eigenvalue over cells =", f"{min(np.linalg.eigvalsh(e).min() for e in effects):.2e}", "(PSD)")
comm = max(np.linalg.norm(effects[i] @ effects[j] - effects[j] @ effects[i], 2)
           for i in range(len(effects)) for j in range(i + 1, len(effects)))
print("   max ||[e_c, e_c']|| =", f"{comm:.4f}", "(compressed code is noncommutative)")
outs = np.zeros(n)
for ci in range(len(uniq)):
    x0 = X[np.where(inv == ci)[0][0]]; hh = x0; A = np.eye(n)
    for W in Ws:
        z = W @ hh; A = np.diag((z > 0).astype(float)) @ W @ A; hh = np.maximum(z, 0)
    outs += A @ effects[ci][1:, 0]
print("   sum_c A_c E[x 1_c] =", np.round(outs, 4), " vs MC E[h_L] =", np.round(h.mean(0), 4))

# ---------------------------------------------------------------------------
# (d) Closure defect accumulation with depth: ARC's mean propagation (K=1,
#     tracking the mean and the trace of the covariance) versus Monte Carlo,
#     at the warm-up shape width 256, depth 16.  Variance-normalised MSE.
# ---------------------------------------------------------------------------
def mean_prop(Ws):
    n = Ws[0].shape[0]; mu = np.zeros(n); s2mean = 1.0; out = []
    for W in Ws:
        mu = W @ mu; s2 = s2mean * np.sum(W ** 2, axis=1); s = np.sqrt(s2)
        a = mu * norm.cdf(mu / s) + s * norm.pdf(mu / s)
        b = (mu ** 2 + s2) * norm.cdf(mu / s) + mu * s * norm.pdf(mu / s)
        out.append(a.copy()); mu = a; s2mean = np.mean(b - a ** 2)
    return out

n, L, N = 256, 16, 1_000_000
Ws = [(rng.normal(size=(n, n)) * np.sqrt(2 / n)).astype(np.float32) for _ in range(L)]
sums = [np.zeros(n) for _ in range(L)]; sq = [np.zeros(n) for _ in range(L)]
chunk = 50_000
for s in range(0, N, chunk):
    hcur = rng.normal(size=(chunk, n)).astype(np.float32)
    for l, W in enumerate(Ws):
        hcur = np.maximum(hcur @ W.T, 0)
        sums[l] += hcur.sum(0, dtype=np.float64); sq[l] += (hcur.astype(np.float64) ** 2).sum(0)
mc = [s_ / N for s_ in sums]; var = [q / N - m ** 2 for q, m in zip(sq, mc)]
mp = mean_prop([W.astype(np.float64) for W in Ws])
print(f"\n(d) width {n}, depth {L}: mean propagation (K=1) vs Monte Carlo ({N} samples)")
print("   layer   MSE(mean-prop)   avg neuron var   MSE/var    MC noise floor (var/N)")
for l in range(L):
    mse = np.mean((mp[l] - mc[l]) ** 2); v = np.mean(var[l])
    print(f"   {l + 1:5d}   {mse:.3e}        {v:.3e}       {mse / v:.2e}   {v / N:.1e}")
