"""TCT-0: zeroth order of the tropical-curvature expansion = Gaussian field + exact one-wall cone lift.

State per layer: pre-activation mean mu (n) and covariance C (n x n) (the thermalised upstream fan, P3 of DESIGN.md).
Lift through the layer's own walls by exact Gaussian cone measures (I3):
  m_j      = s_j phi(h_j) + mu_j Phi(h_j)                     (one wall)
  K_jk     = E[relu(z_j) relu(z_k)]  (two walls: bivariate truncated-normal moments via Owen's T)
Transport: mu' = m W, C' = W^T (K - m m^T) W.
"""
import numpy as np
from scipy.special import owens_t, ndtr
from scipy.stats import norm

SQ2PI = np.sqrt(2 * np.pi)


def phi2(h, k, rho):
    """P(U < h, V < k), standard bivariate normal with correlation rho (vectorised, Owen's T)."""
    h = np.where(np.abs(h) < 1e-12, 1e-12, h); k = np.where(np.abs(k) < 1e-12, 1e-12, k)
    r = np.sqrt(np.maximum(1 - rho * rho, 1e-15))
    t = 0.5 * ndtr(h) + 0.5 * ndtr(k) - owens_t(h, (k - rho * h) / (h * r)) - owens_t(k, (h - rho * k) / (k * r))
    beta = np.where(h * k < 0, 0.5, 0.0)
    return t - beta


def relu_pair(mu1, s1, mu2, s2, rho):
    """E[relu(X) relu(Y)], (X, Y) bivariate normal (means mu, sds s, correlation rho)."""
    h = mu1 / s1; k = mu2 / s2
    r = np.sqrt(np.maximum(1 - rho * rho, 1e-15))
    P = phi2(h, k, rho)
    pa = norm.pdf(h); pb = norm.pdf(k)
    Qa = ndtr((k - rho * h) / r); Qb = ndtr((h - rho * k) / r)
    EU = pa * Qa + rho * pb * Qb
    EV = pb * Qb + rho * pa * Qa
    # E[U V 1(U > -h, V > -k)]
    EUV = rho * P - rho * h * pa * Qa - rho * k * pb * Qb + r / SQ2PI * norm.pdf(np.sqrt(np.maximum((h * h - 2 * rho * h * k + k * k) / (r * r), 0)))
    return s1 * s2 * (EUV + k * EU + h * EV + h * k * P)


def tct0(Ws):
    """Ws: list of (n, n) weights (x @ W convention). Returns (L, n) per-layer means."""
    n = Ws[0].shape[0]
    mu = np.zeros(n); C = Ws[0].T.astype(np.float64) @ Ws[0].astype(np.float64)
    out = []
    for l, W in enumerate(Ws):
        W = W.astype(np.float64)
        if l > 0:
            mu = m @ W; C = W.T @ Cov @ W
        s = np.sqrt(np.maximum(np.diag(C), 1e-300)); h = mu / s
        m = s * norm.pdf(h) + mu * ndtr(h)
        out.append(m)
        if l == len(Ws) - 1:
            break
        rho = np.clip(C / np.outer(s, s), -1 + 1e-12, 1 - 1e-12)
        K = relu_pair(mu[:, None], s[:, None], mu[None, :], s[None, :], rho)
        np.fill_diagonal(K, (mu * mu + s * s) * ndtr(h) + mu * s * norm.pdf(h))
        Cov = K - np.outer(m, m)
    return np.array(out)


if __name__ == "__main__":
    # unit check of relu_pair against Monte Carlo
    rng = np.random.default_rng(0)
    for _ in range(4):
        mu1, mu2 = rng.standard_normal(2); s1, s2 = rng.uniform(.5, 2, 2); rho = rng.uniform(-.95, .95)
        Z = rng.standard_normal((4_000_000, 2)); X = mu1 + s1 * Z[:, 0]; Y = mu2 + s2 * (rho * Z[:, 0] + np.sqrt(1 - rho ** 2) * Z[:, 1])
        print(f"formula {relu_pair(mu1, s1, mu2, s2, rho):.5f}  MC {np.mean(np.maximum(X,0)*np.maximum(Y,0)):.5f}")
