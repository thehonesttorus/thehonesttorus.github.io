"""Stage 21: the gate channel, its cut register, graded mixing, localisation, traces (notes/stage21).

Conventions of this repository: W[l] has shape (n_out, n_in) and z_l = W[l] h_(l-1); the note writes W_l^T for
our W[l]. So the channel of layer l is Phi_l(X) = E[D_l W[l] X W[l]^T D_l] = K_l o (W[l] X W[l]^T), the dual
F_l = Phi_l^*(I) = W[l]^T diag(Phi(alpha_l)) W[l] = sum_j Phi(alpha_j) w_j w_j^T (w_j = rows of W[l]).

  relu_corr_map, gap_closed, gap_quadrature, trajectory, radii     Thm 3.1(a), Cor 3.2
  gate_hermite, cut_kernel, second_moment                           Thm 2.1 (Price channel), Prop 2.2 (level-2 marginal)
  channel, channel_dual, flatness                                   Prop 2.2, Prop 2.3 (trace preservation is flatness)
  level_gain                                                         Thm 3.1(b)
  szegedy_phases                                                     Prop 4.1
  orth_probe_trace, orth_probe_var                                   Thm 4.3
  closure_states                                                     the Gaussian-part trajectory of pass 1
"""
import numpy as np
from scipy.special import ndtr
from whest.relu_gauss import hermite_relu, relu_moments

SQ2PI = np.sqrt(2 * np.pi)


# ---------------------------------------------------------------- Sec 3: graded gap along the arc-cosine trajectory
def relu_corr_map(r):
    """f(rho) = (sqrt(1 - rho^2) + (pi - arccos rho) rho) / pi, the ReLU correlation map."""
    r = np.asarray(r, float); return (np.sqrt(1 - r * r) + (np.pi - np.arccos(r)) * r) / np.pi


def gap_closed(r):
    """g = 1/2 + arcsin(rho)/pi = f'(rho)."""
    return 0.5 + np.arcsin(r) / np.pi


def gap_quadrature(r, deg=200):
    """2 E Phi(alpha)^2 with alpha ~ N(0, rho/(1 - rho)), by Gauss-Hermite quadrature (independent check of (a))."""
    x, w = np.polynomial.hermite_e.hermegauss(deg); sd = np.sqrt(r / (1 - r))
    return 2 * np.sum(w * ndtr(sd * x) ** 2) / np.sqrt(2 * np.pi)


def trajectory(L=16):
    """rho_l (rho_0 = 0) and g_l = f'(rho_(l-1)) for layers l = 1..L."""
    rho = [0.0]
    for _ in range(L): rho.append(float(relu_corr_map(rho[-1])))
    return np.array(rho), np.array([gap_closed(rho[l - 1]) for l in range(1, L + 1)])


def radii(g, k, eta):
    """Cor 3.2: number of most recent layers to keep so that older incoherent level-k information is suppressed
    in amplitude below eta (prod over kept layers of g_l^(k/2)); None if more than len(g)."""
    amp, ell = 1.0, 0
    for gl in g[::-1]:
        amp *= gl ** (k / 2); ell += 1
        if amp < eta: return ell
    return None


# ---------------------------------------------------------------- Sec 2: the Price channel and its cut register
def gate_hermite(a, K):
    """Hermite coefficients of the gate 1(u + a > 0), u ~ N(0, 1): c_0 = Phi(a), c_k = He_(k-1)(-a) phi(a) = d_(k+1)."""
    d = hermite_relu(a, K + 1); return d[1:]


def cut_kernel(mu, C, K=10):
    """Level-2 marginal of the cut law: K_ij = Pr(z_i > 0, z_j > 0) for z ~ N(mu, C) (tetrachoric series), K_ii = Phi(a_i)."""
    sd = np.sqrt(np.clip(np.diag(C), 1e-300, None)); a = mu / sd; c = gate_hermite(a, K)
    rho = C / np.outer(sd, sd); np.fill_diagonal(rho, 0.0)
    Kc = np.outer(c[0], c[0]); rk = np.ones_like(rho); fact = 1.0
    for k in range(1, K + 1):
        rk = rk * rho; fact *= k; Kc += np.outer(c[k], c[k]) * rk / fact
    np.fill_diagonal(Kc, c[0]); return Kc


def second_moment(mu, C, K=10):
    """S^u = E[u u^T] for u = relu(z), z ~ N(mu, C) (the closure's Hermite series)."""
    m, Kh, _, _ = relu_moments(mu, C, K); return Kh + np.outer(m, m), m


def wall_coupling(mu, C):
    """psi_ij = p_(z_i)(0) E[u_j | z_i = 0] (psi_ii = 0): the exact extra response of E[u_i u_j] to a variance change,
    d E[u_i u_j] / d Sigma_ii = psi_ij / 2 for i != j (heat equation: (1/2) E[delta(z_i) u_j]). Thm 2.1's Schur form
    K o dS^z omits it; it carries the mean (to leading order psi_ij = m_j phi(a_i) / sigma_i, cancelling in the centred
    covariance)."""
    from whest.relu_gauss import phi
    d = np.clip(np.diag(C), 1e-300, None); sd = np.sqrt(d)
    cm = mu[None, :] - C * (mu / d)[:, None]                      # E[z_j | z_i = 0]
    cv = np.clip(d[None, :] - C * C / d[:, None], 1e-300, None)    # Var[z_j | z_i = 0]
    cs = np.sqrt(cv); al = cm / cs
    psi = (phi(mu / sd) / sd)[:, None] * (cm * ndtr(al) + cs * phi(al))
    np.fill_diagonal(psi, 0.0); return psi


def full_tangent(mu, C, Kc, dS):
    """Exact tangent of S^u at fixed pre-activation mean: K o dS + (diag(dS) psi + psi^T diag(dS)) / 2."""
    psi = wall_coupling(mu, C); dd = np.diag(dS)
    return Kc * dS + 0.5 * (dd[:, None] * psi + (dd[:, None] * psi).T)


def channel(W, Kc, X):
    """Phi(X) = K o (W X W^T): the tangent of S^u_(l-1) -> S^u_l at fixed means (Thm 2.1)."""
    return Kc * (W @ X @ W.T)


def channel_dual(W, Pa):
    """F = Phi^*(I) = W^T diag(Phi(alpha)) W, so Tr Phi(X) = Tr(F X) (Prop 2.3)."""
    return W.T @ (Pa[:, None] * W)


def flatness(F, mu_prev):
    """Eigenvalues of Pi_perp (F - I) Pi_perp on the complement of the collective direction mu_hat."""
    u = mu_prev / np.linalg.norm(mu_prev); e = np.zeros_like(u); e[0] = 1.0
    v = u - e; v = v / np.linalg.norm(v) if np.linalg.norm(v) > 1e-12 else v
    H = np.eye(len(u)) - 2 * np.outer(v, v); B = H[:, 1:]                 # orthonormal basis of mu_hat^perp
    return np.linalg.eigvalsh(B.T @ (F - np.eye(len(u))) @ B)


def level_gain(W, Pa, V, k):
    """||K^(k) o (W)^(x k) v^(x k)||^2 / ||v||^(2k) for the product (distinct-index) kernel K^(k) = Phi^(x k), i.e.
    (||A W v||^2/||v||^2)^k, for each column v of V."""
    r = np.sum((Pa[:, None] * (W @ V)) ** 2, 0) / np.sum(V * V, 0)
    return r ** k


# ---------------------------------------------------------------- Sec 4: quantum walk, lifting, traces
def szegedy_phases(r, t):
    """Eigenphases +-2 arccos(e^(-r t / 2)) of the reflection walk on Wick level r (Prop 4.1)."""
    return 2 * np.arccos(np.exp(-np.asarray(r) * t / 2))


def orth_probe_trace(A, m, rng):
    """Thm 4.3: (n/m) sum_(i<=m) q_i^T A q_i over the first m columns of a Haar orthogonal matrix."""
    n = A.shape[0]; Q, R = np.linalg.qr(rng.standard_normal((n, n))); Q = Q * np.sign(np.diag(R))[None, :]
    q = Q[:, :m]; return n / m * np.sum(q * (A @ q))


def orth_probe_var(A, m):
    n = A.shape[0]; tr = np.trace(A)
    return 2 * n * (n - m) / (m * (n - 1) * (n + 2)) * (np.sum(A * A) - tr * tr / n)


# ---------------------------------------------------------------- pass 1: the Gaussian-part trajectory
def closure_states(W, K=10, keep_C=True):
    """Forward Schroedinger sweep of the Gaussian part. Per layer l (0-based): pre-activation mean mu, covariance C
    (if keep_C), alpha, Phi(alpha), post-activation mean m and covariance Kh, and rho_l = |m|^2 / Tr S^u."""
    n_in = W[0].shape[1]; mu_h = np.zeros(n_in); C_h = np.eye(n_in); out = []
    for l, Wl in enumerate(W):
        mu = Wl @ mu_h; C = Wl @ C_h @ Wl.T
        m, Kh, sd, a = relu_moments(mu, C, K)
        st = dict(mu=mu, sd=sd, alpha=a, Pa=ndtr(a), m=m, rho=float(m @ m / (np.trace(Kh) + m @ m)),
                  mu_prev=mu_h.copy(), rho_prev=float(mu_h @ mu_h / (np.trace(C_h) + mu_h @ mu_h)))
        if keep_C: st.update(C=C, Kh=Kh, C_prev=C_h)
        out.append(st); mu_h, C_h = m, Kh
    return out
