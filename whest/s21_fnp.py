"""Stage 21 pass 1 written for the scorer (flopscope only, float32, idioms of notes/compute/COST_MODEL.md).

The Gaussian-part forward (Schroedinger) sweep: per layer, the pre-activation mean mu = W m and covariance
C = W K W^T, then the post-activation mean m and covariance K of relu(z), z ~ N(mu, C), by the Hermite series in
the off-diagonal correlation (K_TERMS terms) with the exact diagonal; also the level-1 cut marginal Phi(alpha) and
the realised gap g_l = 2 mean Phi(alpha)^2 (Thm 3.1). Convention: whestbench passes weights W_mlp with
h_next = relu(h @ W_mlp), i.e. our W = W_mlp^T.

Cost-relevant choices (each measured in notes/compute/COST_MODEL.md):
  covariance transport  factor form  L = cholesky(K + eps I), T = W @ L, C = inner(T, T)   3.33 n^3, exactly symmetric
                        (the f32 symmetric-sandwich einsum can raise SymmetryError after being charged)
  Hermite series        Horner in the tagged correlation rho with outer(a_k, a_k) terms        ~14 n^2
  Phi, phi              Abramowitz-Stegun 26.2.17 and exp(-x^2/2), on n-vectors only            O(n)
  no x**k, no where, no float64 scalars; a few dozen counted calls per layer.
"""
import math
import flopscope as flops
import flopscope.numpy as fnp

K_TERMS = 10
SQ2PI_INV = 1.0 / math.sqrt(2.0 * math.pi)
_P = 0.2316419
_B = [c * SQ2PI_INV for c in (0.319381530, -0.356563782, 1.781477937, -1.821255978, 1.330274429)]
_INV_SQRT_FACT = [1.0 / math.sqrt(math.factorial(k)) for k in range(K_TERMS + 2)]


def phi_Phi(x):
    """phi(x), Phi(x) elementwise (abs error 3e-7), sharing one exp."""
    e = fnp.exp(x * x * (-0.5))
    ph = e * SQ2PI_INV
    t = fnp.reciprocal(fnp.abs(x) * _P + 1.0)
    p = t * _B[4] + _B[3]
    for c in (_B[2], _B[1], _B[0]):
        p = p * t + c
    q = e * (p * t)                                    # upper tail Q(|x|)
    return ph, fnp.copysign(0.5 - q, x) + 0.5


def relu_hermite(a, ph, Ph, K):
    """d_0 = phi + a Phi, d_1 = Phi, d_k = phi He_(k-2)(-a) (k >= 2), as a list of n-vectors."""
    x = -a; d = [ph + a * Ph, Ph]
    hm, hc = None, None                                 # He_(j-1), He_j at x
    for k in range(2, K + 1):
        j = k - 2
        if j == 0: hc = x * 0.0 + 1.0
        elif j == 1: hm, hc = hc, x
        else: hm, hc = hc, x * hc - hm * float(j - 1)
        d.append(ph * hc)
    return d


def relu_moments_fnp(mu, C, K=K_TERMS):
    """Mean m and covariance Kh of relu(z), z ~ N(mu, C); C must carry a symmetric tag (an inner() Gram does)."""
    var = fnp.maximum(fnp.diag(C), 1e-30)
    sd = fnp.sqrt(var); a = mu / sd
    ph, Ph = phi_Phi(a)
    d = relu_hermite(a, ph, Ph, K)
    S = fnp.outer(sd, sd)                               # tagged
    rho = C / S
    rho = rho - fnp.diag(fnp.diag(rho))                 # zero diagonal, tag kept
    ak = [d[k] * _INV_SQRT_FACT[k] for k in range(K + 1)]
    acc = fnp.outer(ak[K], ak[K])
    for k in range(K - 1, 0, -1):
        acc = acc * rho + fnp.outer(ak[k], ak[k])
    acc = acc * rho                                     # series starts at k = 1
    Kh = acc * S
    m = sd * d[0]
    second = var * ((a * a + 1.0) * Ph + a * ph)        # E relu(z)^2
    Kh = Kh - fnp.diag(fnp.diag(Kh)) + fnp.diag(second - m * m)
    return m, Kh, a, Ph


def pass1(weights, jitter=1e-6, K=K_TERMS, keep=False):
    """Forward Gaussian-part sweep. weights: list of whestbench (n, n) float32 arrays (h_next = relu(h @ W_mlp)).
    Returns the (L, n) post-activation means, the realised gaps g_l, and (if keep) per-layer states."""
    out, gaps, states = [], [], []
    m = None; Kh = None
    for l, Wm in enumerate(weights):
        W = Wm.T                                        # view: our convention z = W h
        if l == 0:
            mu = fnp.zeros((W.shape[0],), dtype=fnp.float32)
            C = fnp.inner(W, W)                         # input covariance I: C = W W^T, tagged
        else:
            mu = W @ m
            ev = fnp.mean(fnp.diag(Kh)) * jitter
            Lf = fnp.linalg.cholesky(Kh + fnp.diag(fnp.zeros((Kh.shape[0],), dtype=fnp.float32) + ev))
            T = W @ Lf
            C = fnp.inner(T, T)
        m, Kh, a, Ph = relu_moments_fnp(mu, C, K)
        out.append(m); gaps.append(fnp.mean(Ph * Ph) * 2.0)
        if keep: states.append(dict(mu=mu, C=C, m=m, Kh=Kh, alpha=a, Pa=Ph))
    return fnp.stack(out, axis=0), gaps, states
