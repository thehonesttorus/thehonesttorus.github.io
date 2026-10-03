# One-layer checks: transfer identity, Poincare bound with the mean arrow, Theta(1) without it.
import numpy as np
from scipy.stats import norm
rng = np.random.default_rng(0)
n, depth_up, N = 256, 4, 40000
s2 = 2.0 / n                                   # prior variance per weight entry
X = rng.standard_normal((N, n))
H = X
for _ in range(depth_up):                      # upstream state: law of h_{l-1}
    H = np.maximum(H @ (rng.standard_normal((n, n)) * np.sqrt(s2)), 0)
m = H.mean(0); Sig = np.cov(H.T, bias=True)
relu = lambda z: np.maximum(z, 0)
def smooth_relu(a, tau):                       # P_tau relu(a) = E relu(a + sqrt(tau) Z)
    st = np.sqrt(np.maximum(tau, 1e-300)); return a * norm.cdf(a / st) + st * norm.pdf(a / st)
def run(Pbasis, rows=40, draws=300, label=""):
    P = Pbasis @ Pbasis.T if Pbasis.shape[1] else np.zeros((n, n))
    Q = np.eye(n) - P
    Hp = H @ Q                                  # residual activations h'(X)
    tau = s2 * (Hp**2).sum(1)                   # tau(X) = s ||(I-P)h(X)||^2
    var_emp, ann_emp, ann_formula, k1 = [], [], [], []
    for _ in range(rows):
        w = rng.standard_normal(n) * np.sqrt(s2)
        a = H @ (P @ w)
        U = rng.standard_normal((draws, n)) * np.sqrt(s2) @ Q        # residual resamples
        F = relu(a[:, None] + Hp @ U.T).mean(0)                    # F(Pw+u) for each u
        var_emp.append(F.var(ddof=1)); ann_emp.append(F.mean())
        ann_formula.append(smooth_relu(a, tau).mean())             # transfer identity
        psi1 = norm.cdf(a / np.sqrt(tau))                           # (P_tau relu)'(a)
        g = (psi1[:, None] * Hp).mean(0)
        k1.append(s2 * g @ g)                                       # first chaos term
    PSP = Q @ Sig @ Q
    bound = s2 * np.linalg.eigvalsh(PSP).max() / 4
    print(f"[{label}] E Var(f|Pw) (MC over residual) = {np.mean(var_emp):.3e}; first-chaos term = {np.mean(k1):.3e}; "
          f"Poincare bound s2*||P'SigP'||/4 = {bound:.3e}")
    print(f"      transfer identity: max |MC annealed mean - E_X P_tau relu(a)| = {np.max(np.abs(np.array(ann_emp)-np.array(ann_formula))):.2e} "
          f"(MC s.e. ~ {np.sqrt(np.mean(var_emp)/draws):.1e})")
    return np.mean(var_emp)
print(f"n={n}, upstream depth={depth_up}, |m|^2/n={m@m/n:.3f}, ||Sigma||_op={np.linalg.eigvalsh(Sig).max():.3f}, tr Sigma/n={np.trace(Sig)/n:.3f}")
v0 = run(np.zeros((n, 0)), label="V={0}        ")
v1 = run((m / np.linalg.norm(m))[:, None], label="V=span(m)    ")
e1 = np.ones(n) / np.sqrt(n)
v2 = run(e1[:, None], label="V=span(1)    ")
evals, evecs = np.linalg.eigh(Sig)
B = np.linalg.qr(np.column_stack([m, evecs[:, -8:]]))[0]
v3 = run(B, label="V=span(m,top8)")
# target spread across neurons, for scale
W = rng.standard_normal((n, n)) * np.sqrt(s2)
fvals = relu(H[:5000] @ W).mean(0)
print(f"scale: variance of f_j across neurons = {fvals.var():.3e}")
