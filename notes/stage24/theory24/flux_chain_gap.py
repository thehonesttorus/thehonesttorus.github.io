"""The arrow (flux) chain of a single layer, exactly, on small arrangements: how the weights set its gap.

Patterns s in {0,1}^n of z = W x, x ~ N(0, I_d). Stationary law p(s) = P(sign pattern s) (orthant probability of
N(0, C), C = W W^T). Arrows s -> s + e_i carry the stationary flux of the Mehler rotation (Rice):
    c(s, s + e_i) = (1/2pi) P(pattern of z_(-i) = s_(-i) | z_i = 0),   conditional covariance = Schur complement of C.
Generator Q(s, s') = c(s, s')/p(s). Reported: spectral gap of Q, the Cheeger constant bound 1/pi (Gaussian
isoperimetry: c(boundary A) >= (1/sqrt(2 pi)) I(p(A)) >= p(A)/pi for p(A) <= 1/2), the gate influence radius
eta = lambda_max(gate correlation) - 1, and lambda_max of the row-correlation matrix of W, for: orthogonal rows,
square Gaussian W (Marchenko-Pastur rows), a collective (equicorrelated) mode, and near-duplicate rows.
  python notes/stage24/theory24/flux_chain_gap.py"""
import itertools, numpy as np
from scipy.stats import multivariate_normal

rng = np.random.default_rng(7)


def orthant(C, signs):
    """P(signs_k * z_k > 0 for all k), z ~ N(0, C)."""
    if len(signs) == 0: return 1.0
    S = np.diag(signs); Cs = S @ C @ S                                         # P(Sz > 0) = P(-Sz < 0) = cdf of N(0, Cs) at 0
    if len(signs) == 1: return 0.5
    return float(multivariate_normal(mean=np.zeros(len(signs)), cov=Cs, allow_singular=True).cdf(np.zeros(len(signs))))


def chain(C):
    n = len(C); pats = list(itertools.product([0, 1], repeat=n)); idx = {s: k for k, s in enumerate(pats)}
    sg = lambda s: np.array([1.0 if b else -1.0 for b in s])
    p = np.array([orthant(C, sg(s)) for s in pats]); p /= p.sum()
    Q = np.zeros((len(pats), len(pats)))
    for i in range(n):
        o = [k for k in range(n) if k != i]; Sc = C[np.ix_(o, o)] - np.outer(C[o, i], C[i, o]) / C[i, i]
        for rest in itertools.product([0, 1], repeat=n - 1):
            c = orthant(Sc, sg(rest)) / (2 * np.pi)
            s0 = list(rest[:i]) + [0] + list(rest[i:]); s1 = list(rest[:i]) + [1] + list(rest[i:])
            a, b = idx[tuple(s0)], idx[tuple(s1)]; Q[a, b] = c / p[a]; Q[b, a] = c / p[b]
    Q -= np.diag(Q.sum(1))
    Dh = np.sqrt(p); H = (Dh[:, None] * Q) / Dh[None, :]; H = 0.5 * (H + H.T)
    ev = np.sort(np.linalg.eigvalsh(-H)); gap = ev[1]
    # gate correlations
    sgn = np.array([sg(s) for s in pats]); g = (sgn + 1) / 2; m = p @ g; Cg = (g * p[:, None]).T @ g - np.outer(m, m)
    Rg = Cg / np.sqrt(np.outer(np.diag(Cg), np.diag(Cg))); eta = np.linalg.eigvalsh(Rg)[-1] - 1
    R = C / np.sqrt(np.outer(np.diag(C), np.diag(C)))
    return gap, eta, np.linalg.eigvalsh(R)[-1], -Q.diagonal().max(), p.min()


if __name__ == "__main__":
    n = 8
    cases = {"orthogonal rows (independent gates)": np.eye(n)}
    W = rng.standard_normal((n, n)); cases["square Gaussian W (MP rows)"] = W @ W.T
    for beta in (0.3, 0.6, 0.9):
        cases[f"equicorrelated rows beta={beta}"] = (1 - beta) * np.eye(n) + beta * np.ones((n, n))
    for beta in (-0.12, ):
        cases[f"equicorrelated rows beta={beta}"] = (1 - beta) * np.eye(n) + beta * np.ones((n, n))
    W2 = rng.standard_normal((n, n)); W2[1] = W2[0] + 0.05 * rng.standard_normal(n); cases["Gaussian W, rows 0,1 near-duplicate"] = W2 @ W2.T
    print(f"n = {n}: gap of the arrow chain (independent-gate value 2/pi = {2 / np.pi:.4f}; Cheeger h >= 1/pi gives gap >= h^2/(2 q_max))")
    for name, C in cases.items():
        gap, eta, lr, qmax, pmin = chain(C)
        print(f"  {name:40s}: gap {gap:.4f} | gate eta {eta:5.2f} | lambda_max(row corr) {lr:5.2f} | max out-rate {qmax:6.2f} | min p {pmin:.1e}"
              f" | Cheeger lower bound on gap {(1 / np.pi) ** 2 / (2 * qmax):.4f}")
