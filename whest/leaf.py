"""Exact great-circle leaf compiler for a bias-free ReLU network (coupled angular recovery note, Prop. 2 / Thm. 3).

For an orthonormal frame B (n x 2) the circle is x(t) = B (cos t, sin t)^T. On every activation arc each neuron's
pre-activation is a sinusoid a cos t + b sin t, so the whole network is compiled by walking t from 0 to 2 pi and
processing every sign change (event) of every neuron of every layer. An event at (layer l, neuron j) flips one gate;
its effect on all deeper layers is a rank-one update (column vector) x (1 x 2 row), pushed by one matvec per layer.

Returned per leaf: g (n,) = (1/2 pi) * integral over the circle of the last layer's activations (the leaf response of
every output neuron), the number of events per layer, the wall-identity residual (Thm. 3: 2 pi g = sum of derivative
jumps) and whether the gate state closes at t = 2 pi.

Conventions as in whest/relu_gauss.py: W[l] has shape (n_out, n_in), z_l = W[l] @ h_{l-1}.
"""
import numpy as np

TWO_PI = 2.0 * np.pi


def next_roots(Z, t_now, eps=1e-12):
    """For sinusoids a cos t + b sin t (rows of Z), the smallest root in (t_now + eps, 2 pi]; inf if none."""
    a, b = Z[:, 0], Z[:, 1]
    psi = np.arctan2(b, a)                      # a cos t + b sin t = R cos(t - psi), zeros at psi +- pi/2
    r1 = np.mod(psi + np.pi / 2, TWO_PI); r2 = np.mod(psi - np.pi / 2, TWO_PI)
    lo = t_now + eps
    r1 = np.where(r1 > lo, r1, np.inf); r2 = np.where(r2 > lo, r2, np.inf)
    out = np.minimum(r1, r2)
    out[(a == 0) & (b == 0)] = np.inf           # identically zero row: no cut
    return out


def compile_leaf(W, B, max_events=None):
    W = [np.ascontiguousarray(Wl, dtype=np.float64) for Wl in W]
    L = len(W); n = W[0].shape[0]
    B = np.asarray(B, dtype=np.float64)
    # state at t = 0 (cos = 1, sin = 0)
    Z = [None] * L; M = [None] * L; H = B
    for l in range(L):
        Z[l] = W[l] @ H; M[l] = Z[l][:, 0] > 0; H = M[l][:, None] * Z[l]
    M0 = [m.copy() for m in M]
    roots = [next_roots(Z[l], 0.0) for l in range(L)]
    g = np.zeros(n); jumps = np.zeros(n); t_prev = 0.0; counts = np.zeros(L, dtype=np.int64); n_ev = 0
    while True:
        mins = np.array([r.min() for r in roots]); l0 = int(np.argmin(mins)); t_e = mins[l0]
        if not np.isfinite(t_e) or t_e >= TWO_PI:
            break
        j = int(np.argmin(roots[l0]))
        # close the arc [t_prev, t_e] with the current output coefficients
        Hout = M[L - 1][:, None] * Z[L - 1]
        g += Hout[:, 0] * (np.sin(t_e) - np.sin(t_prev)) + Hout[:, 1] * (np.cos(t_prev) - np.cos(t_e))
        dphi_old = -Hout[:, 0] * np.sin(t_e) + Hout[:, 1] * np.cos(t_e)
        # flip the gate: the new value is the sign of the sinusoid just after t_e (its derivative at the root)
        a, b = Z[l0][j]
        new = (-a * np.sin(t_e) + b * np.cos(t_e)) > 0
        sgn = (1.0 if new else 0.0) - (1.0 if M[l0][j] else 0.0)
        M[l0][j] = new
        c = sgn * Z[l0][j]                        # the 1 x 2 coefficient row of the rank-one change of h_{l0}
        v = np.zeros(n); v[j] = 1.0
        for l in range(l0 + 1, L):                # push the rank-one change through the deeper layers
            vz = W[l] @ v
            Z[l] += np.outer(vz, c)
            roots[l] = next_roots(Z[l], t_e)
            v = M[l] * vz
        roots[l0][j] = next_roots(Z[l0][j:j + 1], t_e)[0]
        Hout = M[L - 1][:, None] * Z[L - 1]
        jumps += (-Hout[:, 0] * np.sin(t_e) + Hout[:, 1] * np.cos(t_e)) - dphi_old
        t_prev = t_e; counts[l0] += 1; n_ev += 1
        if max_events is not None and n_ev >= max_events:
            break
    Hout = M[L - 1][:, None] * Z[L - 1]
    g += Hout[:, 0] * (np.sin(TWO_PI) - np.sin(t_prev)) + Hout[:, 1] * (np.cos(t_prev) - np.cos(TWO_PI))
    g /= TWO_PI
    closed = all(np.array_equal(M[l], M0[l]) for l in range(L))
    return dict(g=g, events=counts, wall_residual=float(np.max(np.abs(jumps / TWO_PI - g))), closed=closed)


def random_frame(n, rng):
    Q, _ = np.linalg.qr(rng.standard_normal((n, 2)))
    return Q


def c_radial(n):
    """E|X| for X ~ N(0, I_n): sqrt(2) Gamma((n+1)/2) / Gamma(n/2)."""
    from scipy.special import gammaln
    return np.sqrt(2.0) * np.exp(gammaln((n + 1) / 2) - gammaln(n / 2))
