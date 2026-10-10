"""Test: the arrow chain of a Gaussian sign pattern has its spectral gap on Walsh level 1 (linear functions of gates).

z = mu + W x, x ~ N(0, I), stationary under the Mehler rotation; gates g_i = 1{z_i > 0}. Rice fluxes (one way):
    c_i(rest) = (1/2pi) exp(-alpha_i^2/2) P(sign pattern of z_(-i) = rest | z_i = 0),  alpha_i = mu_i / sd_i,
with the conditional law N(mu_(-i) - C_(-i,i) mu_i / C_ii, Schur complement). Dirichlet form of a linear f = v.g:
    E(f, f) = sum_i v_i^2 Gamma_i,   Gamma_i = sum_rest c_i(rest) = exp(-alpha_i^2/2) / (2 pi),
so the best linear Rayleigh quotient is lambda_lin = 1 / lambda_max(Gamma^(-1/2) Cov(g) Gamma^(-1/2)) >= gap.
The conjecture is gap = lambda_lin. Exact enumeration (n <= 10), Genz orthants at tight tolerance.
  python notes/stage24/theory24/gap_level1.py"""
import itertools, numpy as np
from scipy.stats import multivariate_normal
from scipy.special import ndtr

rng = np.random.default_rng(11)


def orthant(m, C, signs):
    """P(signs_k * z_k > 0 for all k), z ~ N(m, C)."""
    k = len(signs)
    if k == 0: return 1.0
    if k == 1: return float(ndtr(signs[0] * m[0] / np.sqrt(C[0, 0])))
    S = np.diag(signs)                                                    # P(S z > 0) = P(-S z < 0), -S z ~ N(-S m, S C S)
    return float(multivariate_normal(mean=-S @ m, cov=S @ C @ S, allow_singular=True, abseps=1e-9, releps=1e-6, maxpts=200_000 * k).cdf(np.zeros(k)))


def analyse(mu, C):
    n = len(C); pats = list(itertools.product([0, 1], repeat=n)); idx = {s: k for k, s in enumerate(pats)}
    sg = lambda s: np.array([1.0 if b else -1.0 for b in s])
    p = np.array([orthant(mu, C, sg(s)) for s in pats]); tot = p.sum(); p /= tot
    Q = np.zeros((len(pats), len(pats))); sd = np.sqrt(np.diag(C)); al = mu / sd
    for i in range(n):
        o = [k for k in range(n) if k != i]
        Sc = C[np.ix_(o, o)] - np.outer(C[o, i], C[i, o]) / C[i, i]; mc = mu[o] - C[o, i] * mu[i] / C[i, i]
        for rest in itertools.product([0, 1], repeat=n - 1):
            c = np.exp(-0.5 * al[i] ** 2) / (2 * np.pi) * orthant(mc, Sc, sg(rest))
            s0 = rest[:i] + (0,) + rest[i:]; s1 = rest[:i] + (1,) + rest[i:]
            a, b = idx[s0], idx[s1]; Q[a, b] = c / p[a]; Q[b, a] = c / p[b]
    Q -= np.diag(Q.sum(1))
    Dh = np.sqrt(p); H = (Dh[:, None] * Q) / Dh[None, :]; H = 0.5 * (H + H.T); ev, U = np.linalg.eigh(-H); gap = ev[1]
    g = np.array([[1.0 if b else 0.0 for b in s] for s in pats]); m = p @ g; Cg = (g * p[:, None]).T @ g - np.outer(m, m)
    Gam = np.exp(-0.5 * al ** 2) / (2 * np.pi); M = Cg / np.sqrt(np.outer(Gam, Gam)); lin = 1 / np.linalg.eigvalsh(M)[-1]
    # Walsh level of the gap eigenvector: share of its variance explained by linear functions of the gates
    f = U[:, 1] / Dh; f -= p @ f; X = g - m; beta = np.linalg.lstsq(X * Dh[:, None], f * Dh, rcond=None)[0]
    share = np.sum((Dh * (X @ beta)) ** 2) / np.sum((Dh * f) ** 2)
    return gap, lin, share, ev[2], abs(tot - 1)


if __name__ == "__main__":
    cases = []
    for n in (5, 7):
        W = rng.standard_normal((n, n)); cases.append((f"n={n} Gaussian W, mu=0", np.zeros(n), W @ W.T))
        cases.append((f"n={n} Gaussian W, random mu (alpha ~ N(0,1))", rng.standard_normal(n) * np.sqrt(np.diag(W @ W.T)), W @ W.T))
        B = rng.standard_normal((n, 2)); C = 0.3 * np.eye(n) + B @ B.T
        cases.append((f"n={n} two collective modes + noise, mu=0", np.zeros(n), C))
        cases.append((f"n={n} two collective modes, alpha ~ N(0, 2^2)", 2 * rng.standard_normal(n) * np.sqrt(np.diag(C)), C))
        A = rng.standard_normal((n, n)); C = A @ A.T + 0.05 * np.eye(n)
        cases.append((f"n={n} ill-conditioned Wishart, mild means", 0.5 * rng.standard_normal(n) * np.sqrt(np.diag(C)), C))
    n = 8; W = rng.standard_normal((n, n)); cases.append((f"n=8 Gaussian W, random mu", rng.standard_normal(n) * np.sqrt(np.diag(W @ W.T)), W @ W.T))
    print("case | gap | linear bound 1/lambda_max(Gamma^-1/2 Cov(g) Gamma^-1/2) | gap/linear | level-1 share of the gap eigenvector | next eigenvalue | orthant sum error")
    for name, mu, C in cases:
        gap, lin, share, e2, err = analyse(mu, C)
        print(f"  {name:48s} | {gap:.6f} | {lin:.6f} | {gap / lin:.6f} | {share:.6f} | {e2:.4f} | {err:.1e}", flush=True)
