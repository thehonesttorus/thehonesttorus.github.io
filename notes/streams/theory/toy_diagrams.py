"""Generic diagram enumerator for the joint cumulant kappa3(f(z_0), f(z_1), f(z_2)) of a near-Gaussian triple, and the
toy validations of the vertex-weight rule and of every closure2 coefficient.

Expansion (formal, asymptotic): E[F(z)] = E_G[ exp( sum_{m>=3} kappa_m . d^m / m! ) F ] with E_G Gaussian of the same
mean and covariance, and E_G[prod_v f(z_v)] = exp( sum_{u<v} C_uv d_u d_v ) prod_v E[f(z_v)] (derivatives in the means).
Expanding both exponentials, a term is a multiset of hyperedges; a hyperedge is a multiset of vertices of size m >= 2
(m = 2 with two distinct vertices: a C edge; m >= 3: a cumulant; a size-2 multiset on one vertex is the variance, which
stays in the Gaussian marginal).  The joint cumulant is the sum over CONNECTED multisets touching all three vertices of

    prod_h kappa_h / prod_v r_{h,v}!  x  prod_types 1/mult!  x  prod_v w_v(deg v),   w_v(d) = E_G[f^(d)(z_v)].

Nothing here is relu-specific: the enumerator takes the cumulant tensors and the vertex weights as inputs.

    python toy_diagrams.py            all checks (identity, completeness, lambda-scaling, relu, dressing)
"""
import itertools
import os
import sys
from math import factorial, sqrt, pi, erf

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

V = 3


# ---------------------------------------------------------------------------------------------------------------
# enumeration
# ---------------------------------------------------------------------------------------------------------------
def all_types(mmax):
    ts = []
    for m in range(2, mmax + 1):
        for c in itertools.product(range(m + 1), repeat=V):
            if sum(c) != m:
                continue
            if m == 2 and max(c) == 2:
                continue  # variance: in the Gaussian marginal
            ts.append(c)
    return ts


def order_generic(c):
    m = sum(c)
    return 1 if m == 2 else m - 2


def order_network(c):
    """fixed random MLP at width n, eps = n^-1/2: kappa_m ~ eps^(m-2) if every multiplicity is even, else eps^(m-1)
    (C edges: eps)."""
    m = sum(c)
    return m - 2 if all(x % 2 == 0 for x in c) else m - 1


def enumerate_diagrams(N, order=order_generic, mmax=None, allow_self=True):
    """all connected multisets of hyperedge types spanning the 3 vertices with total order <= N.
    returns list of (diagram, order), diagram = tuple of (type, mult)."""
    if mmax is None:
        mmax = N + 2
    types = [c for c in all_types(mmax) if order(c) <= N and order(c) >= 1]
    if not allow_self:
        types = [c for c in types if sum(1 for x in c if x) >= 2]
    out = []

    def rec(i, budget, cur):
        if i == len(types):
            if cur and connected(cur):
                out.append((tuple(cur), N - budget))
            return
        rec(i + 1, budget, cur)
        o = order(types[i])
        mult = 1
        while o * mult <= budget:
            rec(i + 1, budget - o * mult, cur + [(types[i], mult)])
            mult += 1
    rec(0, N, [])
    return out


def degrees(diag):
    d = [0] * V
    for c, m in diag:
        for v in range(V):
            d[v] += c[v] * m
    return d


def connected(diag):
    if min(degrees(diag)) < 1:
        return False
    parent = list(range(V))

    def find(x):
        while parent[x] != x:
            x = parent[x]
        return x
    for c, m in diag:
        vs = [v for v in range(V) if c[v]]
        for a in vs[1:]:
            parent[find(a)] = find(vs[0])
    return len({find(v) for v in range(V)}) == 1


def canon(diag):
    best = None
    for p in itertools.permutations(range(V)):
        d2 = []
        for c, m in diag:
            c2 = [0] * V
            for v in range(V):
                c2[p[v]] = c[v]
            d2.append((tuple(c2), m))
        d2 = tuple(sorted(d2))
        if best is None or d2 < best:
            best = d2
    return best


def weight(diag, kap, w):
    """kap: dict sorted-index-tuple -> cumulant (incl. pairs for C); w[v][d] vertex weights."""
    val = 1.0
    for c, m in diag:
        idx = tuple(sorted(v for v in range(V) for _ in range(c[v])))
        h = kap[idx] / np.prod([factorial(x) for x in c])
        val *= h ** m / factorial(m)
    for v, d in enumerate(degrees(diag)):
        val *= w[v][d]
    return val


# ---------------------------------------------------------------------------------------------------------------
# toy model: z = mu + L g + t Q(g), g ~ N(0, I_3)
# ---------------------------------------------------------------------------------------------------------------
class Toy:
    def __init__(self, s, t, seed=0, mu=(0.3, -0.2, 0.1), sig=(1.0, 0.9, 1.1), cubic=0.0):
        rng = np.random.default_rng(seed)
        A = rng.standard_normal((V, V)); np.fill_diagonal(A, 0.0)
        self.L = np.diag(sig) @ (np.eye(V) + s * A)
        B = rng.standard_normal((V, V, V)); B = 0.5 * (B + B.transpose(0, 2, 1))
        self.B = 0.4 * B
        G3 = rng.standard_normal((V, V)) * 0.3
        self.mu = np.array(mu); self.t = t; self.cubic = cubic; self.G3 = G3

    def z(self, g):
        """g: (P, V) -> z (P, V)."""
        q = np.einsum("pa,vab,pb->pv", g, self.B, g) - np.trace(self.B, axis1=1, axis2=2)[None, :]
        z = self.mu + g @ self.L.T + self.t * q
        if self.cubic:
            z = z + self.cubic * self.t ** 2 * ((g ** 3 - 3 * g) @ self.G3.T)
        return z


def gh_grid(q):
    x, w = np.polynomial.hermite_e.hermegauss(q)
    w = w / w.sum()
    X = np.stack(np.meshgrid(x, x, x, indexing="ij"), -1).reshape(-1, V)
    Wq = (w[:, None, None] * w[None, :, None] * w[None, None, :]).ravel()
    return X, Wq


def cumulants(toy, mmax, q=16):
    """exact joint cumulants of z (all sorted index tuples of length 2..mmax) by Gauss-Hermite quadrature
    (exact: the integrands are polynomials in g of degree <= 2 mmax (6 mmax with the cubic term) < 2q)."""
    X, Wq = gh_grid(q)
    Z = toy.z(X)
    U = Z - Wq @ Z
    mom = {}

    def cm(idx):
        if idx not in mom:
            mom[idx] = float(Wq @ np.prod(U[:, list(idx)], axis=1)) if idx else 1.0
        return mom[idx]
    kap = {}
    for m in range(2, mmax + 1):
        for idx in itertools.combinations_with_replacement(range(V), m):
            pos = list(range(m))
            tot = 0.0
            for part in partitions_no_singletons(pos):
                b = len(part)
                tot += (-1) ** (b - 1) * factorial(b - 1) * np.prod([cm(tuple(sorted(idx[p] for p in blk))) for blk in part])
            kap[idx] = tot
    mean = Wq @ Z
    return kap, mean


def partitions_no_singletons(pos):
    if not pos:
        yield []
        return
    first, rest = pos[0], pos[1:]
    for r in range(1, len(rest) + 1):
        for comb in itertools.combinations(rest, r):
            remaining = [p for p in rest if p not in comb]
            for p in partitions_no_singletons(remaining):
                yield [[first] + list(comb)] + p


def f_smooth(z, delta):
    """Gaussian-smoothed relu E_xi[relu(z + delta xi)] (entire); E_N(mu,s2)[f^(d)] = relu weights at variance s2 + delta^2."""
    if delta == 0:
        return np.maximum(z, 0.0)
    u = z / delta
    return z * 0.5 * (1 + erf_v(u / sqrt(2))) + delta * np.exp(-0.5 * u * u) / sqrt(2 * pi)


def erf_v(x):
    try:
        from scipy.special import erf as e
        return e(x)
    except ImportError:
        return np.vectorize(erf)(x)


def kappa3_f(toy, delta, q):
    X, Wq = gh_grid(q)
    F = f_smooth(toy.z(X), delta)
    Fc = F - Wq @ F
    return float(Wq @ (Fc[:, 0] * Fc[:, 1] * Fc[:, 2])), Wq @ (toy.z(X) > 0).astype(float)


def relu_weights(mu, var, dmax):
    from numpy.polynomial.hermite_e import hermeval
    sig = sqrt(var); al = mu / sig
    phi = np.exp(-0.5 * al ** 2) / sqrt(2 * pi)
    Phi = 0.5 * (1 + erf(al / sqrt(2)))
    w = [sig * (al * Phi + phi), Phi]
    for d in range(2, dmax + 1):
        w.append(float(hermeval(-al, [0] * (d - 2) + [1]) * phi / sig ** (d - 1)))
    return w


def vertex_weights(kap, mean, delta, dmax=40):
    return [relu_weights(mean[v], kap[(v, v)] + delta ** 2, dmax) for v in range(V)]


# ---------------------------------------------------------------------------------------------------------------
# closure2 terms as diagram shapes (vertex 0 = i, 1 = j, 2 = k)
# ---------------------------------------------------------------------------------------------------------------
E_ij, E_ik, E_jk = (1, 1, 0), (1, 0, 1), (0, 1, 1)
SHAPES = {
    "WICK": [[(E_ik, 1), (E_jk, 1)]],
    "B0": [[((1, 1, 1), 1)]],
    "B1": [[((2, 0, 1), 1), (E_jk, 1)]],
    "B2": [[((2, 0, 1), 1), (E_ij, 1)]],
    "B3": [[((3, 0, 0), 1), (E_ij, 1), (E_ik, 1)]],
    "B4": "C3",
    "B5": [[((2, 2, 0), 1), (E_ik, 1)]],
    "B6": [[((2, 1, 1), 1)]],
    "B7": [[((1, 1, 1), 1), (E_ij, 1)]],
    "S1": "C4",
    "S2c": [[((4, 0, 0), 1), (E_ij, 1), (E_ik, 1)]],
    "S3a1": [[((1, 1, 1), 1), (E_ij, 2)]],
    "S3a2": [[((1, 1, 1), 1), (E_ij, 1), (E_ik, 1)]],
    "S3b1": [[((2, 0, 1), 1), (E_jk, 2)]],
    "S3b2": [[((2, 0, 1), 1), (E_ij, 2)]],
    "S3b3": [[((2, 0, 1), 1), (E_ij, 1), (E_jk, 1)]],
    "S3b4": [[((2, 0, 1), 1), (E_ij, 1), (E_ik, 1)]],
    "S3b5": [[((2, 0, 1), 1), (E_jk, 1), (E_ik, 1)]],
    "S4a": [[((1, 1, 1), 2)]],
    "S4b": [[((1, 1, 1), 1), ((2, 0, 1), 1)]],
    "S4c": [[((2, 0, 1), 1), ((0, 2, 1), 1)]],
    "S4d": [[((1, 0, 2), 1), ((0, 1, 2), 1)]],
    "S4e": [[((2, 0, 1), 1), ((0, 1, 2), 1)]],
    "S5a": [[((2, 1, 1), 1), (E_jk, 1)]],
    "S5b": [[((2, 1, 1), 1), (E_ij, 1)]],
    "S5c": [[((3, 1, 0), 1), (E_ik, 1)]],
    "S5d": [[((3, 1, 0), 1), (E_jk, 1)]],
    "S6a": [[((2, 2, 0), 1), (E_ik, 2)]],
    "S6b": [[((2, 2, 0), 1), (E_ik, 1), (E_jk, 1)]],
    "S6c": [[((2, 2, 0), 1), (E_ik, 1), (E_ij, 1)]],
    "S7a": [[((2, 2, 0), 1), ((1, 1, 1), 1)]],
    "S7b": [[((2, 2, 0), 1), ((2, 0, 1), 1)]],
    "S7c": [[((2, 2, 0), 1), ((1, 0, 2), 1)]],
    "S8": [[((2, 2, 0), 1), ((0, 2, 2), 1)]],
    "S9a": [[((2, 2, 1), 1)]],
    "S9b": [[((3, 1, 1), 1)]],
    "S10": [[((2, 2, 2), 1)]],
    "G3t": [[(E_ij, 1), (E_jk, 1), (E_ik, 1)]],
    "G4a": [[(E_ij, 2), (E_jk, 2)]],
    "G4b": [[(E_ij, 2), (E_jk, 1), (E_ik, 1)]],
}
NONLEAF = ["B0", "B6", "B7", "G3t", "S3a1", "S3a2", "S3b1", "S3b2", "S3b3", "S4a", "S4b", "S4c", "S4d", "S4e", "S5a",
           "S5b", "S6a", "S6b", "S7a", "S7b", "S7c", "S8", "G4a", "G4b", "S9a", "S9b", "S10"]


def has_c_leaf(diag):
    """some vertex has non-self degree 1 and that leg is on a C edge."""
    nonself = [(c, m) for c, m in diag if sum(1 for x in c if x) >= 2]
    dn = degrees(nonself)
    for v in range(V):
        if dn[v] == 1:
            for c, m in nonself:
                if c[v] and sum(c) == 2:
                    return True
    return False


def check_leaf_completeness(N=4):
    """network order <= N: every diagram is a C-leaf diagram (resummed by Phi_k C_ik Gamma_ij) or, up to gate
    dressings of degree-1 vertices, one of the NONLEAF shapes; and no NONLEAF shape has a C-leaf."""
    diags = enumerate_diagrams(N, order_network, mmax=N + 2)
    nl = set()
    for name in NONLEAF:
        nl |= shape_signatures(name, [d for d, _ in diags])
    bad_nl = [nm for nm in NONLEAF if any(has_c_leaf(sg) for sg in shape_signatures(nm, [d for d, _ in diags]))]
    missing = []
    nleaf = 0
    for d, o in diags:
        if has_c_leaf(d):
            nleaf += 1
            continue
        nonself = tuple(sorted((c, m) for c, m in d if sum(1 for x in c if x) >= 2))
        if canon(nonself) in nl:
            continue
        missing.append((o, canon(d)))
    print(f"leaf decomposition at network order <= {N}: {nleaf} C-leaf diagrams, {len(diags) - nleaf - len(missing)} non-leaf "
          f"covered by {len(NONLEAF)} shapes, uncovered: {len(missing)}; NONLEAF shapes with a C-leaf: {bad_nl or 'none'}")
    for o, sg in missing[:10]:
        print("   ", o, sg)
    return missing


def shape_signatures(name, diags):
    sh = SHAPES[name]
    if isinstance(sh, str):
        ne = int(sh[1:])
        return {canon(d) for d in diags if all(sum(c) == 2 for c, m in d) and sum(m for c, m in d) == ne}
    return {canon(tuple(sorted(d))) for d in sh}


def toy_objects_n3(kap, mean, delta=0.0):
    """closure2's object dict for the toy triple (n = 3), with BARE (Gaussian) degree-1 weights."""
    import closure2 as c2
    n = V
    C = np.zeros((n, n)); T = np.zeros((n, n, n)); K4f = np.zeros((n, n, n))
    P = np.zeros((n, n, n)); Q = np.zeros((n, n, n)); S = np.zeros((n, n, n))
    for i, j in itertools.product(range(n), repeat=2):
        C[i, j] = kap[tuple(sorted((i, j)))]
    for i, j, k in itertools.product(range(n), repeat=3):
        T[i, j, k] = kap[tuple(sorted((i, j, k)))]
        K4f[i, j, k] = kap[tuple(sorted((i, i, j, k)))]
        if (i, j, k) in [(a, b, c) for a, b, c in itertools.permutations(range(3))]:
            P[i, j, k] = kap[tuple(sorted((i, i, j, j, k)))]
            Q[i, j, k] = kap[tuple(sorted((i, i, i, j, k)))]
            S[i, j, k] = kap[tuple(sorted((i, i, j, j, k, k)))]
    var = np.diag(C).copy()
    w = c2.relu_w(mean, var + delta ** 2)
    return dict(mu=mean, var=var, C=C, T=T, K4f=K4f, P=P, Q=Q, S=S, Phi=w[1], w=w)


# ---------------------------------------------------------------------------------------------------------------
# checks
# ---------------------------------------------------------------------------------------------------------------
def check_closure2_identity(s=0.35, t=0.3, verbose=True):
    """every closure2 term (coef x tensor at the all-distinct entry) against the enumerator's sum over all labelled
    diagrams of the same shape, on a non-Gaussian toy triple with relu vertex weights (an algebraic identity)."""
    import closure2 as c2
    toy = Toy(s, t, seed=1, cubic=1.0)
    kap, mean = cumulants(toy, 6, q=24)
    o = toy_objects_n3(kap, mean)
    diags = [d for d, _ in enumerate_diagrams(4, order_generic, mmax=6)]
    w = [relu_weights(mean[v], kap[(v, v)], 16) for v in range(V)]
    # gather enumerator sums per signature
    by_sig = {}
    for d in diags:
        by_sig.setdefault(canon(d), 0.0)
        by_sig[canon(d)] += weight(d, kap, w)
    TT = c2.terms(o)
    worst = 0.0
    lines = []
    for name in [k for k in SHAPES if k != "WICK"]:
        if name not in TT:
            continue
        sigs = shape_signatures(name, diags)
        ref = sum(by_sig.get(sg, 0.0) for sg in sigs)
        val = c2.TERM_INFO[name][2] * TT[name][0, 1, 2]
        err = abs(val - ref) / max(abs(ref), 1e-300)
        worst = max(worst, err)
        lines.append(f"  {name:5s} coef {c2.TERM_INFO[name][2]:6.4f}  closure2 {val:+.10e}  enumerator {ref:+.10e}  rel.diff {err:.1e}")
    # the oracle's B3 coefficient 1.0 against the enumerator
    sigB3 = shape_signatures("B3", diags)
    refB3 = sum(by_sig.get(sg, 0.0) for sg in sigB3)
    lines.append(f"  B3 with oracle_k3.CLOSURE_COEF = 1.0: {1.0 * TT['B3'][0, 1, 2]:+.6e} vs enumerator {refB3:+.6e} (ratio {TT['B3'][0,1,2]/refB3:.4f})")
    # Wick leading term
    import oracle_k3 as ok
    Kw = ok.wick_model(o["C"], o["Phi"], o["w"][2], np.zeros((3, 3, 3)))
    refW = sum(by_sig.get(sg, 0.0) for sg in shape_signatures("WICK", diags))
    lines.append(f"  WICK  closure2/oracle {Kw[0,1,2]:+.10e}  enumerator {refW:+.10e}  rel.diff {abs(Kw[0,1,2]-refW)/abs(refW):.1e}")
    if verbose:
        print(f"closure2 term identity on a toy triple (s={s}, t={t}): worst relative difference {worst:.2e}")
        print("\n".join(lines))
    return worst


def check_completeness(N=4):
    """network-order enumeration: every connected diagram of order <= N is either a closure2 term, part of the Wick /
    hermite terms, or a self-hyperedge decoration of a degree-1 vertex (exact through the true gate P(z > 0))."""
    diags = enumerate_diagrams(N, order_network, mmax=N + 2)
    known = {}
    for name in SHAPES:
        for sg in shape_signatures(name, [d for d, _ in diags]):
            known[sg] = name
    missing = {}
    gate_dressed = 0
    for d, o in diags:
        sg = canon(d)
        if sg in known:
            continue
        # strip self-hyperedges sitting on vertices whose non-self degree is 1
        nonself = [(c, m) for c, m in d if sum(1 for x in c if x) >= 2]
        dn = degrees(nonself) if nonself else [0] * V
        selfs = [(c, m) for c, m in d if sum(1 for x in c if x) == 1]
        ok_ = all(dn[[v for v in range(V) if c[v]][0]] == 1 for c, m in selfs)
        if ok_ and selfs and nonself and canon(tuple(sorted(nonself))) in known:
            gate_dressed += 1
            continue
        missing.setdefault(o, []).append(sg)
    print(f"completeness at network order <= {N}: {len(diags)} connected diagrams; {gate_dressed} are degree-1 gate dressings "
          f"(exact through P(z>0)); not covered by a closure2 term: " + (", ".join(f"order {k}: {len(v)}" for k, v in sorted(missing.items())) or "none"))
    for k, v in sorted(missing.items()):
        for sg in v[:12]:
            print(f"    order {k}: {sg}")
    return missing


def check_scaling(delta=0.5, lams=(0.4, 0.2, 0.1, 0.05, 0.025), N=5, q=64, s0=1.0, t0=1.0, cubic=1.0):
    """exact kappa3(f(z)) of the toy (quadrature) against the generic-order truncations of the enumerator, with
    s = s0 lam, t = t0 lam: the residual after order N must scale as lam^(N+1)."""
    diags = enumerate_diagrams(N, order_generic, mmax=N + 2)
    print(f"\nlambda-scaling, smoothed relu delta={delta}: {len(diags)} diagrams up to generic order {N}; "
          f"rows lam, exact kappa3, |residual| after order 1..{N} (relative to exact)")
    res = []
    for lam in lams:
        toy = Toy(s0 * lam, t0 * lam, seed=2, cubic=cubic)
        kap, mean = cumulants(toy, N + 2, q=4 * (N + 2))
        w = vertex_weights(kap, mean, delta, dmax=2 * N + 8)
        ex, _ = kappa3_f(toy, delta, q)
        part = np.zeros(N + 1)
        for d, o in diags:
            part[o] += weight(d, kap, w)
        cum = np.cumsum(part)
        r = [abs(ex - cum[k]) for k in range(1, N + 1)]
        res.append(r)
        print(f"  lam {lam:6.3f}  exact {ex:+.6e}  " + " ".join(f"{x / abs(ex):.2e}" for x in r))
    res = np.array(res)
    slopes = np.log(res[:-1] / res[1:]) / np.log(np.array(lams[:-1]) / np.array(lams[1:]))[:, None]
    print("  local slopes d log|res| / d log lam (expected N+1 = " + ", ".join(str(k + 1) for k in range(1, N + 1)) + "):")
    for i in range(len(lams) - 1):
        print(f"    lam {lams[i]:.3f}->{lams[i+1]:.3f}: " + " ".join(f"{x:5.2f}" for x in slopes[i]))
    return res, slopes


def check_closure_ladder_relu(lams=(0.3, 0.2, 0.1), q=(160, 224), N=5):
    """relu itself (no smoothing): exact kappa3 by high-order quadrature (two grids; their difference bounds the
    accuracy), against the enumerator's generic-order truncations with relu vertex weights."""
    diags = enumerate_diagrams(N, order_generic, mmax=N + 2)
    print(f"\nrelu toy triple (delta = 0), generic-order truncations 1..{N} vs quadrature on {q[0]}^3 and {q[1]}^3 grids")
    for lam in lams:
        toy = Toy(lam, lam, seed=2, cubic=1.0)
        kap, mean = cumulants(toy, N + 2, q=4 * (N + 2))
        exs = [kappa3_f(toy, 0.0, qq)[0] for qq in q]
        ex = exs[-1]; acc = abs(exs[-1] - exs[0])
        w = vertex_weights(kap, mean, 0.0, dmax=2 * N + 8)
        part = np.zeros(N + 1)
        for d, o in diags:
            part[o] += weight(d, kap, w)
        cum = np.cumsum(part)
        print(f"  lam {lam:.2f}: exact {ex:+.6e} (grid diff {acc / abs(ex):.1e} rel.)  rel. residual after order 1..{N}: "
              + " ".join(f"{abs(ex - cum[k]) / abs(ex):.2e}" for k in range(1, N + 1)))


def check_dressing(lams=(0.4, 0.2, 0.1, 0.05, 0.025), delta=0.5, q=64):
    """the dressed vertex weight E_true[f^(d)(z_v)] (d = 1, 2) of the smoothed relu, exact by quadrature, against the
    Gaussian weight plus the self-hyperedge series kappa3/3! w(d+3) + kappa4/4! w(d+4) + kappa3^2/(2 3!^2) w(d+6)
    (orders t^0, t^1, t^2).  At delta -> 0 these are P(z > 0) and the density of z at 0."""
    print(f"\ndressing of a vertex (smoothed relu delta={delta}): |E f^(d)(z) - series| after orders 0, 1, 2; d = 1 and d = 2")
    X, Wq = gh_grid(q)
    for lam in lams:
        toy = Toy(lam, lam, seed=3, cubic=1.0)
        kap, mean = cumulants(toy, 4, q=16)
        z0 = toy.z(X)[:, 0]
        u = z0 / delta
        e1 = float(Wq @ (0.5 * (1 + erf_v(u / sqrt(2)))))                 # E f'(z)
        e2 = float(Wq @ (np.exp(-0.5 * u * u) / sqrt(2 * pi) / delta))     # E f''(z)
        w = relu_weights(mean[0], kap[(0, 0)] + delta ** 2, 12)
        k3, k4 = kap[(0, 0, 0)], kap[(0, 0, 0, 0)]
        row = []
        for d, ex in ((1, e1), (2, e2)):
            s0 = w[d]; s1 = s0 + k3 / 6 * w[d + 3]; s2 = s1 + k4 / 24 * w[d + 4] + k3 ** 2 / 72 * w[d + 6]
            row.append(f"d={d}: {abs(ex - s0):.2e} {abs(ex - s1):.2e} {abs(ex - s2):.2e}")
        print(f"  lam {lam:.3f}  " + "   ".join(row))


def check_scaling_table(delta=0.5, lams=(0.2, 0.1, 0.07, 0.05, 0.035, 0.025), N=5, q=72):
    """signed residual after generic order N divided by lam^(N+1): must tend to a constant as lam -> 0."""
    diags = enumerate_diagrams(N, order_generic, mmax=N + 2)
    print(f"\nresidual_N / lam^(N+1), smoothed relu delta={delta}, quadrature {q}^3 ({len(diags)} diagrams up to order {N})")
    print("   lam    " + " ".join(f"   N={k}    " for k in range(1, N + 1)))
    for lam in lams:
        toy = Toy(lam, lam, seed=2, cubic=1.0)
        kap, mean = cumulants(toy, N + 2, q=4 * (N + 2))
        w = vertex_weights(kap, mean, delta, dmax=2 * N + 8)
        ex, _ = kappa3_f(toy, delta, q)
        part = np.zeros(N + 1)
        for d, o in diags:
            part[o] += weight(d, kap, w)
        cum = np.cumsum(part)
        print(f"  {lam:6.3f}  " + " ".join(f"{(ex - cum[k]) / lam ** (k + 1):+.4e}" for k in range(1, N + 1)))


def check_leaf_identity(delta=0.5, lams=(0.2, 0.1, 0.07, 0.05, 0.035, 0.025), N=5, q=72):
    """exact kappa3 - [sum_roles Phi_k C_ik Gamma_ij - sum_centers Phi_j Phi_k C_ij C_ik p_i] (all from quadrature:
    Gamma_ij = Cov(f'(z_i), f(z_j)), Phi = E f', p = E f'') must equal the sum of the diagrams WITHOUT a C-leaf;
    the residual after the non-leaf diagrams of generic order <= N, divided by lam^(N+1), must tend to a constant."""
    diags = [(d, o) for d, o in enumerate_diagrams(N, order_generic, mmax=N + 2) if not has_c_leaf(d)]
    X, Wq = gh_grid(q)
    print(f"\nleaf identity, smoothed relu delta={delta}: (exact - leaf + two-leaf - non-leaf diagrams to order N) / lam^(N+1)")
    print("   lam     |leaf part|/|k3|  " + " ".join(f"   N={k}    " for k in range(1, N + 1)))
    for lam in lams:
        toy = Toy(lam, lam, seed=2, cubic=1.0)
        kap, mean = cumulants(toy, N + 2, q=4 * (N + 2))
        w = vertex_weights(kap, mean, delta, dmax=2 * N + 8)
        Z = toy.z(X)
        F = f_smooth(Z, delta); F1 = 0.5 * (1 + erf_v(Z / delta / sqrt(2))); F2 = np.exp(-0.5 * (Z / delta) ** 2) / sqrt(2 * pi) / delta
        Fc = F - Wq @ F
        ex = float(Wq @ (Fc[:, 0] * Fc[:, 1] * Fc[:, 2]))
        Phi = Wq @ F1; p = Wq @ F2
        F1c = F1 - Phi
        Gam = np.einsum("p,pi,pj->ij", Wq, F1c, Fc)
        C = np.array([[kap[tuple(sorted((i, j)))] for j in range(V)] for i in range(V)])
        leaf = 0.0
        for k, i, j in itertools.permutations(range(V)):
            leaf += Phi[k] * C[i, k] * Gam[i, j]
        two = 0.0
        for i in range(V):
            j, k = [v for v in range(V) if v != i]
            two += Phi[j] * Phi[k] * C[i, j] * C[i, k] * p[i]
        part = np.zeros(N + 1)
        for d, o in diags:
            part[o] += weight(d, kap, w)
        cum = np.cumsum(part)
        rem = ex - leaf + two
        print(f"  {lam:6.3f}   {abs(leaf - two) / abs(ex):9.3f}       " + " ".join(f"{(rem - cum[k]) / lam ** (k + 1):+.4e}" for k in range(1, N + 1)))


def in_tclass(diag):
    """exactly one all-distinct kappa3 hyperedge, every other non-self hyperedge on one and the same pair."""
    nonself = [(c, m) for c, m in diag if sum(1 for x in c if x) >= 2]
    tri = [(c, m) for c, m in nonself if all(c)]
    if len(tri) != 1 or tri[0][0] != (1, 1, 1) or tri[0][1] != 1:
        return False
    pairs = {tuple(v for v in range(V) if c[v]) for c, m in nonself if not all(c)}
    return len(pairs) <= 1


def check_tclass_identity(delta=0.5, lams=(0.2, 0.1, 0.07, 0.05, 0.035, 0.025), N=5, q=72):
    """exact - LEAF + TWOLEAF - T_ijk (Phi_i Phi_j Phi_k + sum_pairs Phi_k Cov(f'(z_i), f'(z_j))) must equal the sum of
    the diagrams that are neither C-leaf nor in the T class; residual after order N over lam^(N+1)."""
    diags = [(d, o) for d, o in enumerate_diagrams(N, order_generic, mmax=N + 2) if not has_c_leaf(d) and not in_tclass(d)]
    X, Wq = gh_grid(q)
    print(f"\nT-class + leaf identity, smoothed relu delta={delta}: residual after the remaining diagrams to order N, / lam^(N+1)")
    for lam in lams:
        toy = Toy(lam, lam, seed=2, cubic=1.0)
        kap, mean = cumulants(toy, N + 2, q=4 * (N + 2))
        w = vertex_weights(kap, mean, delta, dmax=2 * N + 8)
        Z = toy.z(X)
        F = f_smooth(Z, delta); F1 = 0.5 * (1 + erf_v(Z / delta / sqrt(2))); F2 = np.exp(-0.5 * (Z / delta) ** 2) / sqrt(2 * pi) / delta
        Fc = F - Wq @ F
        ex = float(Wq @ (Fc[:, 0] * Fc[:, 1] * Fc[:, 2]))
        Phi = Wq @ F1; p = Wq @ F2
        F1c = F1 - Phi
        Gam = np.einsum("p,pi,pj->ij", Wq, F1c, Fc)
        Gg = np.einsum("p,pi,pj->ij", Wq, F1c, F1c)
        C = np.array([[kap[tuple(sorted((i, j)))] for j in range(V)] for i in range(V)])
        leaf = sum(Phi[k] * C[i, k] * Gam[i, j] for k, i, j in itertools.permutations(range(V)))
        two = sum(Phi[[v for v in range(V) if v != i][0]] * Phi[[v for v in range(V) if v != i][1]]
                  * C[i, [v for v in range(V) if v != i][0]] * C[i, [v for v in range(V) if v != i][1]] * p[i] for i in range(V))
        T = kap[(0, 1, 2)]
        tcl = T * (Phi[0] * Phi[1] * Phi[2] + Phi[2] * Gg[0, 1] + Phi[0] * Gg[1, 2] + Phi[1] * Gg[0, 2])
        part = np.zeros(N + 1)
        for d, o in diags:
            part[o] += weight(d, kap, w)
        cum = np.cumsum(part)
        rem = ex - leaf + two - tcl
        print(f"  {lam:6.3f}  |T-class|/|k3| {abs(tcl) / abs(ex):6.3f}   " + " ".join(f"{(rem - cum[k]) / lam ** (k + 1):+.4e}" for k in range(1, N + 1)))


if __name__ == "__main__":
    check_closure2_identity()
    check_completeness(4)
    check_leaf_completeness(4)
    check_leaf_identity()
    check_tclass_identity()
    check_dressing()
    check_scaling_table(delta=0.5)
    check_scaling_table(delta=0.25, q=112)
