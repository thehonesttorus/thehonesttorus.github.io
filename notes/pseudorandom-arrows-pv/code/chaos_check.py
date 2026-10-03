# Exact check of the chaos formula  Var(F|Pw) = sum_k (1/k!) E_{X,X'}[psi_k psi_k' K^k]
import numpy as np
from scipy.stats import norm
from scipy.special import eval_hermitenorm, factorial
rng = np.random.default_rng(1)
n, Ns = 256, 3000
s2 = 2.0 / n
H = rng.standard_normal((Ns, n))
for _ in range(4):
    H = np.maximum(H @ (rng.standard_normal((n, n)) * np.sqrt(s2)), 0)
m = H.mean(0); e = m / np.linalg.norm(m)
Q = np.eye(n) - np.outer(e, e)
Hp = H @ Q; tau = s2 * (Hp**2).sum(1); st = np.sqrt(tau)
K = s2 * Hp @ Hp.T                                  # K(X,X') on the empirical law (incl. diagonal)
for trial in range(3):
    w = rng.standard_normal(n) * np.sqrt(s2)
    a = H @ (np.outer(e, e) @ w)
    # psi_k = (P_tau relu)^{(k)}(a):  k=1: Phi(a/st); k>=2: (-1)^{k} He_{k-2}(a/st) phi(a/st) / st^{k-1}
    terms = []
    for k in range(1, 13):
        if k == 1: psi = norm.cdf(a / st)
        else: psi = (-1)**k * eval_hermitenorm(k - 2, a / st) * norm.pdf(a / st) / st**(k - 1)
        terms.append((psi @ (K**k) @ psi) / Ns**2 / factorial(k))
    U = rng.standard_normal((20000, n)) * np.sqrt(s2) @ Q
    F = np.maximum(a[:, None] + Hp @ U.T, 0).mean(0)
    print(f"row {trial}: chaos sum (k<=12) = {sum(terms):.4e}  [k=1: {terms[0]:.2e}, k=2: {terms[1]:.2e}, k=3: {terms[2]:.2e}, k=4: {terms[3]:.2e}]"
          f"   MC Var over 20000 residual draws = {F.var(ddof=1):.4e} (+/- {F.var(ddof=1)*np.sqrt(2/20000):.1e})")
