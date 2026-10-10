"""Euler-Stein defect of a chain output with the first-jet (frozen-covariance) Laplacian.

Every hidden unit of a bias-free ReLU net is positively homogeneous of degree 1 in the input, so its mean under
N(t, I) satisfies  m(t) = Delta_t m(t)  (Euler + Stein).  For a chain output e_l (all layers), the first-jet Laplacian
is the exact Laplacian of the mean-only chain with the covariance frozen at the untilted run:
    J_0 = Phi_0 (.) W_1,  Delta m_0 = m_0,
    A_l = W_l J_{l-1},  J_l = Phi_l (.) A_l,  Delta m_l = (phi_l / sigma_l) (.) rowsumsq(A_l) + Phi_l (.) (W_l Delta m_{l-1}).
Fields r_l, sigma_l are recovered from the chain's own output (m = sigma G(r), mu = W m_{l-1}).
Reports, per layer: corr(defect, error), best merge constant a, MSE before/after (a = best, 1/2, 1/3, 1/4).

    python scripts/es_jlap.py DATA NET OUT_NPY [--selftest]
"""
import sys, numpy as np
from scipy.stats import norm

G = lambda r: norm.pdf(r) + r * norm.cdf(r)
rg = np.linspace(-9, 9, 36001); qq = rg / G(rg)                       # q = mu/m = r/G(r), increasing in r


def fields(m, mu):
    r = np.interp(mu / np.maximum(m, 1e-300), qq, rg)
    s = m / np.maximum(G(r), 1e-300)
    return r, s


def jlap(W, out, W0_exact=True):
    """First-jet Laplacian of every layer's mean, from the chain's own per-layer means `out` (L, n)."""
    L, n = out.shape
    lap = np.zeros_like(out); Phi = np.zeros_like(out); R = np.zeros_like(out); S = np.zeros_like(out)
    w0 = np.linalg.norm(W[0], axis=1)
    R[0] = 0.0; S[0] = w0; Phi[0] = 0.5
    J = 0.5 * W[0]; lap[0] = norm.pdf(0) * w0                                # = m_0 exactly
    for l in range(1, L):
        mu = W[l] @ out[l - 1]
        r, s = fields(out[l], mu); R[l], S[l] = r, s; Phi[l] = norm.cdf(r)
        A = W[l] @ J
        lap[l] = (norm.pdf(r) / s) * np.einsum("ij,ij->i", A, A) + Phi[l] * (W[l] @ lap[l - 1])
        J = Phi[l][:, None] * A
    return lap, R, S, Phi


def mean_chain(W, t, S, Phi0=None):
    """Mean-only chain with frozen sigmas S: m_l(t)."""
    L, n = S.shape; out = np.zeros_like(S)
    mu = W[0] @ t; out[0] = S[0] * G(mu / S[0])
    for l in range(1, L):
        mu = W[l] @ out[l - 1]; out[l] = S[l] * G(mu / S[l])
    return out


def selftest():
    rng = np.random.default_rng(0); n, L = 48, 4
    W = rng.standard_normal((L, n, n)) * np.sqrt(2 / n)
    S = np.ones((L, n)); S[0] = np.linalg.norm(W[0], axis=1)
    for l in range(1, L): S[l] = 0.8 + 0.3 * rng.random(n)
    out = mean_chain(W, np.zeros(n), S)
    lap, _, _, _ = jlap(W, out)
    h = 1e-3; num = np.zeros_like(out)
    for i in range(n):
        e = np.zeros(n); e[i] = h
        num += (mean_chain(W, e, S) + mean_chain(W, -e, S) - 2 * out) / h ** 2
    print(f"selftest: max |symbolic - probe Laplacian| = {np.abs(lap - num).max():.2e} (scale {np.abs(lap).max():.2f})")


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        selftest(); sys.exit()
    D, net, outp = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    W = np.load(f"{D}/W_off{net}.npy").astype(np.float64); truth = np.load(f"{D}/truth_off{net}.npz")["m"]
    out = np.load(outp).astype(np.float64)
    lap, R, S, Phi = jlap(W, out)
    e = out - truth; d = out - lap
    L = out.shape[0]
    print(f"net {net}: {outp}")
    print(" l   MSE        |defect| rms  corr(d,e)  a*     MSE(a*)    MSE(1/2)   MSE(1/3)   MSE(1/4)   |lap|rms  <r> ")
    for l in range(L):
        el, dl = e[l], d[l]
        c = np.corrcoef(dl, el)[0, 1]; a = (dl @ el) / (dl @ dl)
        mse = lambda aa: np.mean((el - aa * dl) ** 2)
        print(f"{l:2d} {np.mean(el**2):.3e}  {np.sqrt(np.mean(dl**2)):.3e}    {c:+.3f}   {a:+.3f}  {mse(a):.3e}  {mse(.5):.3e}  {mse(1/3):.3e}  {mse(.25):.3e}  {np.sqrt(np.mean(lap[l]**2)):.3f}  {R[l].mean():+.2f}")
    # the layerwise defect transported to the output by the first jet, as a second regressor
    T = lambda l, v: Phi[l] * (W[l] @ v)
    cols = []
    for l in range(L):
        v = d[l].copy()
        for k in range(l + 1, L): v = T(k, v)
        cols.append(v)
    X = np.array(cols).T; y = e[-1]
    coef, *_ = np.linalg.lstsq(X, y, rcond=None); res = y - X @ coef
    print(f"last layer: regress error on all transported layer defects: R^2 = {1 - res @ res / (y @ y):.3f}; coefs {np.round(coef, 3)}")
