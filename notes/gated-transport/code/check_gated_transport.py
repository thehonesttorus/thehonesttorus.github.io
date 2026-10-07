# Exact checks of the gated transport at the independent reference (note XXXIV).   python check_gated_transport.py
# Reference: z_a ~ N(mu_a, s_a^2) independent, y_a = relu(z_a). A first-order tangent of the law of z on index support S with
# multiplicities k_a (a pre-activation cumulant entry) has Hermite score  omega * prod_a He_(k_a)  in precision coordinates,
# score_1 = q, score_2 = q^2 - 1/s^2, score_3 = q^3 - 3 q/s^2,  q = (z - mu)/s^2.
# The check computes the first-order change of a post-activation joint cumulant DIRECTLY: raw moments of y under p(1 + eps h)
# by 1-D quadrature of score x power-of-relu (factorised over independent coordinates), cumulants by the partition formula,
# exact linearisation in eps. No Stein identity is used. It is compared with the closed forms of the note:
#   c(p, k) = E[(d/dz)^k (y - m)^p]  (distributional),  and for the classes the chain drops:
#   (3,1):   d kappa(y_a,y_a,y_a,y_b) = omega [c_a(3,k) - 3 v_a c_a(1,k)] c_b(1,l)
#   (2,1,1): d kappa(y_a,y_a,y_b,y_c) = omega c_a(2,k_a) c_b(1,k_b) c_c(1,k_c)
#   (2,2):   d kappa(y_a,y_a,y_b,y_b) = omega c_a(2,k) c_b(2,l)
#   (1,1,1,1): d kappa(y_a,y_b,y_c,y_d) = omega prod c(1,k)
# with omega = C_ab (covariance), Gamma_ab/2 (kappa3 pair entry), T_abc (all-distinct kappa3), K_ab/4 ((2,2) kappa4 entry),
# B_ab/6 ((3,1) kappa4 entry), V/2 ((2,1,1) kappa4 entry): the multinomial count of the entry's orderings over 4! or 3!.
import itertools, math, numpy as np
from numpy.polynomial.legendre import leggauss
X, Wt = leggauss(400)
def E(f, mu, s):
    lo, hi = mu - 12 * s, mu + 12 * s
    cuts = [lo] + ([0.0] if lo < 0 < hi else []) + [hi]
    tot = 0.0
    for a, b in zip(cuts[:-1], cuts[1:]):
        z = 0.5 * (b - a) * X + 0.5 * (b + a)
        tot += 0.5 * (b - a) * np.sum(Wt * f(z) * np.exp(-0.5 * ((z - mu) / s) ** 2) / (s * math.sqrt(2 * math.pi)))
    return tot
def score(k, mu, s):
    q = lambda z: (z - mu) / s ** 2
    return {0: lambda z: np.ones_like(z), 1: q, 2: lambda z: q(z) ** 2 - 1 / s ** 2,
            3: lambda z: q(z) ** 3 - 3 * q(z) / s ** 2}[k]
relu = lambda z: np.maximum(z, 0.0)

def dcumulant(idx, prm, tang):
    """first-order change of kappa(y_idx[0], ..., y_idx[-1]) under the tangent tang = (omega, {neuron: k})."""
    omega, ks = tang
    def raw(block, lin):                                  # E[prod_{i in block} y_i] (lin=False) or its eps-coefficient (lin=True)
        pw = {}
        for i in block:
            pw[i] = pw.get(i, 0) + 1
        neurons = set(pw) | (set(ks) if lin else set())
        out = omega if lin else 1.0
        for a in neurons:
            mu, s = prm[a]; p = pw.get(a, 0); k = ks.get(a, 0) if lin else 0
            out *= E(lambda z: relu(z) ** p * score(k, mu, s)(z), mu, s)
        return out
    def partitions(seq):
        if not seq:
            yield []
            return
        first, rest = seq[0], seq[1:]
        for p in partitions(rest):
            yield [[first]] + p
            for i in range(len(p)):
                yield p[:i] + [[first] + p[i]] + p[i + 1:]
    pos = list(range(len(idx)))
    tot = 0.0
    for part in partitions(pos):
        coef = (-1) ** (len(part) - 1) * math.factorial(len(part) - 1)
        blocks = [[idx[i] for i in B] for B in part]
        base = [raw(b, False) for b in blocks]; lin = [raw(b, True) for b in blocks]
        tot += coef * sum(lin[j] * np.prod([base[t] for t in range(len(blocks)) if t != j]) for j in range(len(blocks)))
    return tot

def c(p, k, mu, s):                                       # E[(d/dz)^k (y - m)^p], closed forms
    a = mu / s; Ph = 0.5 * (1 + math.erf(a / math.sqrt(2))); ph = math.exp(-a * a / 2) / math.sqrt(2 * math.pi)
    m = mu * Ph + s * ph; Ey2 = (mu * mu + s * s) * Ph + mu * s * ph; v = Ey2 - m * m
    dl, dl1 = ph / s, -a * ph / s ** 2                   # E delta(z), E delta'(z)
    tab = {(1, 1): Ph, (1, 2): dl, (1, 3): dl1,
           (2, 1): 2 * m * (1 - Ph), (2, 2): 2 * Ph - 2 * m * dl,
           (3, 1): 3 * (v - m * m * (1 - Ph)), (3, 2): 6 * m * (1 - Ph) + 3 * m * m * dl,
           (3, 3): 6 * Ph - 6 * m * dl + 3 * m * m * dl1}
    return tab[(p, k)], v

prm = {0: (0.7, 1.3), 1: (-0.4, 0.8), 2: (0.2, 1.1), 3: (-1.0, 0.9)}
cases = []
# (3,1)_y on {0,1}: kappa(y0,y0,y0,y1)
for nm, om, ks in [("C_01", 0.13, {0: 1, 1: 1}), ("Gamma_01", 0.21 / 2, {0: 2, 1: 1}), ("Gamma_10", -0.17 / 2, {0: 1, 1: 2}),
                   ("K_01", 0.09 / 4, {0: 2, 1: 2}), ("B_01", 0.11 / 6, {0: 3, 1: 1}), ("B_10", -0.05 / 6, {0: 1, 1: 3})]:
    k, l = ks[0], ks[1]
    ca, va = c(3, k, *prm[0]); c1a, _ = c(1, k, *prm[0]); cb, _ = c(1, l, *prm[1])
    cases.append((f"(3,1) from {nm}", dcumulant([0, 0, 0, 1], prm, (om, ks)), om * (ca - 3 * va * c1a) * cb))
# (2,2)_y from a kappa3 pair entry and from a (2,2) kappa4 entry
for nm, om, ks in [("Gamma_01", 0.21 / 2, {0: 2, 1: 1}), ("K_01", 0.09 / 4, {0: 2, 1: 2}), ("C_01", 0.13, {0: 1, 1: 1})]:
    pred = om * c(2, ks[0], *prm[0])[0] * c(2, ks[1], *prm[1])[0]
    cases.append((f"(2,2) from {nm}", dcumulant([0, 0, 1, 1], prm, (om, ks)), pred))
# (2,1,1)_y on {0,1,2}: kappa(y0,y0,y1,y2)
for nm, om, ks in [("T_012", 0.19, {0: 1, 1: 1, 2: 1}), ("V_(0;12)", 0.07 / 2, {0: 2, 1: 1, 2: 1}), ("V_(1;02)", 0.07 / 2, {0: 1, 1: 2, 2: 1})]:
    pred = om * c(2, ks[0], *prm[0])[0] * c(1, ks[1], *prm[1])[0] * c(1, ks[2], *prm[2])[0]
    cases.append((f"(2,1,1) from {nm}", dcumulant([0, 0, 1, 2], prm, (om, ks)), pred))
# (1,1,1,1)_y on {0,1,2,3}
om = 0.05
cases.append(("(1,1,1,1) from U_0123", dcumulant([0, 1, 2, 3], prm, (om, {0: 1, 1: 1, 2: 1, 3: 1})),
              om * np.prod([c(1, 1, *prm[a])[0] for a in range(4)])))
# support preservation: unary and partial-support tangents leave the mixed cumulants unchanged
cases.append(("(3,1) from unary kappa3 tau_0", dcumulant([0, 0, 0, 1], prm, (0.3 / 6, {0: 3})), 0.0))
cases.append(("(2,1,1) from pair Gamma_01", dcumulant([0, 0, 1, 2], prm, (0.21 / 2, {0: 2, 1: 1})), 0.0))
cases.append(("(2,1,1) from C_01", dcumulant([0, 0, 1, 2], prm, (0.13, {0: 1, 1: 1})), 0.0))
for nm, got, pred in cases:
    print(f"{nm:28s} quadrature {got:+.10f}  closed form {pred:+.10f}  diff {abs(got - pred):.1e}")
# the cumulant theory's N(0, I) witness: T_112 = t gives t(3c/8 + 3c^3/2) and t(3c/8 - 3c^3/2)
p0 = {0: (0.0, 1.0), 1: (0.0, 1.0)}; cc = 1 / math.sqrt(2 * math.pi)
print("N(0,I) witness:", dcumulant([0, 0, 0, 1], p0, (0.5, {0: 2, 1: 1})), "vs", 3 * cc / 8 + 1.5 * cc ** 3,
      "|", dcumulant([1, 1, 1, 0], p0, (0.5, {0: 2, 1: 1})), "vs", 3 * cc / 8 - 1.5 * cc ** 3)

# --- Ward consistency: the radial tangent of a Gaussian state, pushed through the gated transport, is the radial tangent of y.
# At z ~ N(mu, diag s^2) the mean-preserving radial tangent (note XXXII section 4, per unit v) has entries
#   C_ab = mu_a mu_b,  Gamma_ab = kappa3-tangent(a,a,b) = 2 mu_b s_a^2,  K_ab = 4 s_a^2 s_b^2  (and unary parts),
# and the exact image is the radial tangent of y: d Cov(y_a,y_b) = m_a m_b, d kappa(y_a,y_a,y_b,y_b) = 4 v_a v_b,
# d kappa(y_a,y_a,y_a,y_b) = 3 m_b kappa3(y_a), d kappa(y_a,y_a,y_b,y_c) = 0 (y's covariance and kappa3 are diagonal here).
print("\nWard consistency of the gated transport (radial tangent in, radial tangent out):")
def radial_tangent_terms(a, b):
    mua, sa = prm[a]; mub, sb = prm[b]
    return [(mua * mub, {a: 1, b: 1}), (2 * mub * sa ** 2 / 2, {a: 2, b: 1}), (2 * mua * sb ** 2 / 2, {a: 1, b: 2}),
            (4 * sa ** 2 * sb ** 2 / 4, {a: 2, b: 2})]
def ystats(a):
    mu, s = prm[a]
    m = E(lambda z: relu(z), mu, s); v = E(lambda z: (relu(z) - m) ** 2, mu, s); k3 = E(lambda z: (relu(z) - m) ** 3, mu, s)
    return m, v, k3
for idx, name, exact in [([0, 1], "Cov(y0,y1)", lambda: ystats(0)[0] * ystats(1)[0]),
                         ([0, 0, 1, 1], "kappa(y0,y0,y1,y1)", lambda: 4 * ystats(0)[1] * ystats(1)[1]),
                         ([0, 0, 0, 1], "kappa(y0,y0,y0,y1)", lambda: 3 * ystats(1)[0] * ystats(0)[2]),
                         ([1, 1, 1, 0], "kappa(y1,y1,y1,y0)", lambda: 3 * ystats(0)[0] * ystats(1)[2])]:
    got = sum(dcumulant(idx, prm, t) for t in radial_tangent_terms(0, 1))
    print(f"  {name:20s} gated transport {got:+.10f}   radial tangent of y {exact():+.10f}   diff {abs(got - exact()):.1e}")
