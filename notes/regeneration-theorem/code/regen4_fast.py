# Note XLIX: one-step regeneration of the fourth-order pair slices of z' = W relu(z) from the true layer-l state,
# by automatic enumeration of the linked-cluster structures (Lemma 1), kept by coherence power counting (Lemma 2).
#   python regen4.py NET L [MCDIR] [WDIR]      transition L -> L+1
# Site factors: shift derivatives of the site cumulants of relu under the true (Edgeworth) marginal of z_a.
# C-edges: the true covariance (off-diagonal). Hyperedges: the true cross-site slices D21, K22, K31 of layer L.
# Targets: the layer-(L+1) truth (full + halves for noise-free R^2) on the active set, readout-weighted.
import sys, re, itertools, math, numpy as np
from scipy.special import ndtr

net, L = int(sys.argv[1]), int(sys.argv[2])
MC = sys.argv[3] if len(sys.argv) > 3 else "."
WD = sys.argv[4] if len(sys.argv) > 4 else "."
EMAX = 4
import os
SLICES = os.environ.get("SLICES", "31,22,diag").split(",")
GSF = os.environ.get("GSF", "0") == "1"
EXPO_MIN = int(os.environ.get("EXPO_MIN", "1"))
EXACT_FALLBACK = os.environ.get("EXACT_FALLBACK", "0") == "1"   # 1: evaluate cyclic coincidence terms by einsum (n^4)
NPAIRS = int(os.environ.get("NPAIRS", "20000"))   # Lemma 2 cut (1 = leading only; 0 adds the n^(-1/2) structures)   # also evaluate with Gaussian site factors (consistency check)
f64 = lambda a: np.asarray(a, dtype=np.float64)
T = {h: np.load(f"{MC}/mc2_off{net}_{h}.npz") for h in ("full", "h0", "h1")}
W = f64(np.load(f"{WD}/W_off{net}.npy", mmap_mode="r")[L + 1])
n = W.shape[0]
off = ~np.eye(n, dtype=bool)
SQ2PI = math.sqrt(2 * math.pi)


def he(k, x):
    h0, h1 = np.ones_like(x), x
    if k == 0:
        return h0
    for m in range(1, k):
        h0, h1 = h1, x * h1 - m * h0
    return h1


def site_factors(mu, var, k3, k4, edgeworth=True, D=4):
    """F[(p, d)] = d-th shift derivative of the p-th cumulant of relu(z_a), p = 1..4, d = 0..D."""
    s = np.sqrt(var); r = mu / s
    g3 = k3 / s**3 if edgeworth else 0 * mu
    g4 = k4 / s**4 if edgeworth else 0 * mu
    a = {0: np.ones_like(mu), 3: g3 / 6, 4: g4 / 24, 6: g3 * g3 / 72}
    ph = np.exp(-0.5 * r * r) / SQ2PI
    x = -r
    Phi = ndtr(r) + ph * sum(c * he(k - 1, x) for k, c in a.items() if k >= 1)
    # p_t^(m)(x) = (-1)^m phi(x) sum a_k He_(k+m)(x); p_z^(m)(0) = p_t^(m)(-r) / s^(m+1)
    pz = [((-1) ** m) * ph * sum(c * he(k + m, x) for k, c in a.items()) / s ** (m + 1) for m in range(4)]
    dR0 = [Phi, pz[0], -pz[1], pz[2], -pz[3]]          # shift derivatives of P(z > 0)
    # raw moments R_p = E relu^p by quadrature over t in [-r, -r + span]
    G = 6001
    span = np.maximum(12.0, 12.0 - x)
    u = np.linspace(0.0, 1.0, G)[None, :]
    t = x[:, None] + span[:, None] * u
    dens = np.exp(-0.5 * t * t) / SQ2PI * sum(c[:, None] * he(k, t) for k, c in a.items())
    z = mu[:, None] + s[:, None] * t
    wq = np.full(G, 2.0); wq[1::2] = 4.0; wq[0] = wq[-1] = 1.0
    wq = wq / (3.0 * (G - 1))
    R = {0: Phi}
    for p in range(1, 5):
        R[p] = (z**p * dens) @ wq * span
    # Taylor polynomials in the shift delta, truncated at order D
    def dR(p, d):
        if d <= p:
            return math.perm(p, d) * R[p - d]
        return math.factorial(p) * dR0[d - p]
    poly = {p: np.stack([dR(p, d) / math.factorial(d) for d in range(D + 1)], axis=1) for p in range(1, 5)}
    def mul(A, B):
        out = np.zeros_like(A)
        for i in range(D + 1):
            out[:, i:] += A[:, i:i + 1] * B[:, :D + 1 - i]
        return out
    R1, R2, R3, R4 = poly[1], poly[2], poly[3], poly[4]
    R1s = mul(R1, R1)
    k = {1: R1, 2: R2 - R1s,
         3: R3 - 3 * mul(R2, R1) + 2 * mul(R1s, R1),
         4: R4 - 4 * mul(R3, R1) - 3 * mul(R2, R2) + 12 * mul(R2, R1s) - 6 * mul(R1s, R1s)}
    return {(p, d): k[p][:, d] * math.factorial(d) for p in range(1, 5) for d in range(D + 1)}


def set_partitions(items):
    items = list(items)
    if not items:
        yield []
        return
    first, rest = items[0], items[1:]
    for part in set_partitions(rest):
        yield [[first]] + part
        for i in range(len(part)):
            yield part[:i] + [[first] + part[i]] + part[i + 1:]


def compositions(total, k):
    if k == 1:
        yield (total,)
        return
    for x in range(total + 1):
        for rest in compositions(total - x, k - 1):
            yield (x,) + rest


def edge_multisets(slots, emax):
    """All multisets of C-edges over vertex-pair slots with total <= emax: dict slot -> multiplicity."""
    def rec(i, left):
        if i == len(slots):
            yield {}
            return
        for m in range(left + 1):
            for rest in rec(i + 1, left - m):
                d = dict(rest)
                if m:
                    d[slots[i]] = m
                yield d
    yield from rec(0, emax)


def connected(nv, links):
    par = list(range(nv))
    def f(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    for grp in links:
        for u in grp[1:]:
            par[f(u)] = f(grp[0])
    return len({f(v) for v in range(nv)}) == 1


def exponent(legrows, verts, vsite, nsite, emult_site):
    """Lemma 2: f - e, maximised over pairings (C-only structures)."""
    p = [0] * nsite; q = [0] * nsite
    for v, legs in enumerate(verts):
        for lg in legs:
            if legrows[lg] == "i":
                p[vsite[v]] += 1
            else:
                q[vsite[v]] += 1
    ident = [q[s] > 0 or p[s] % 2 == 1 for s in range(nsite)]
    changed = True
    while changed:
        changed = False
        for (s1, s2), m in emult_site.items():
            if m % 2 == 1 and not (ident[s1] and ident[s2]):
                ident[s1] = ident[s2] = True; changed = True
    nid = sum(ident); nfree = nsite - nid
    e = sum(emult_site.values())
    return nid + 2 * nfree - e


def exponent_h(legrows, verts, vsite, nsite, emult_site, hsites, coherent):
    """Lemma 2 with one hyperedge on sites hsites: random-signed (identifies both) or coherent (self-pairable)."""
    p = [0] * nsite; q = [0] * nsite
    for v, legs in enumerate(verts):
        for lg in legs:
            if legrows[lg] == "i":
                p[vsite[v]] += 1
            else:
                q[vsite[v]] += 1
    ident = [q[s] > 0 or p[s] % 2 == 1 for s in range(nsite)]
    if not coherent:
        ident[hsites[0]] = ident[hsites[1]] = True
    changed = True
    while changed:
        changed = False
        for (s1, s2), m in emult_site.items():
            if m % 2 == 1 and not (ident[s1] and ident[s2]):
                ident[s1] = ident[s2] = True; changed = True
    nid = sum(ident); nfree = nsite - nid
    return nid + 2 * nfree - sum(emult_site.values())


def structures(legrows, hyper):
    """Yield (verts, vsite, nsite, cedges{(v,w):m}, hedge or None, weight, expo)."""
    legs = list(range(len(legrows)))
    for P in set_partitions(legs):
        nv = len(P)
        for S in set_partitions(range(nv)):
            vsite = [0] * nv
            for si, blk in enumerate(S):
                for v in blk:
                    vsite[v] = si
            ns = len(S)
            slots = [(v, w) for v in range(nv) for w in range(v + 1, nv) if vsite[v] != vsite[w]]
            # C-only structures
            for E in edge_multisets(slots, EMAX):
                if not connected(nv, [list(k) for k in E]):
                    continue
                es = {}
                for (v, w), m in E.items():
                    key = tuple(sorted((vsite[v], vsite[w])))
                    es[key] = es.get(key, 0) + m
                expo = exponent(legrows, P, vsite, ns, es)
                wgt = 1.0
                for m in E.values():
                    wgt /= math.factorial(m)
                yield P, vsite, ns, E, None, wgt, expo
            # one carried hyperedge + at most one C-edge
            if not hyper or ns < 2:
                continue
            for name, (r1, r2) in (("D21", (2, 1)), ("K31", (3, 1)), ("K22", (2, 2))):
                for X in range(ns):
                    for Y in range(ns):
                        if X == Y or (name == "K22" and Y < X):
                            continue
                        vX = [v for v in range(nv) if vsite[v] == X]
                        vY = [v for v in range(nv) if vsite[v] == Y]
                        for cx in compositions(r1, len(vX)):
                            for cy in compositions(r2, len(vY)):
                                inc = {}
                                for v, m in zip(vX, cx):
                                    if m:
                                        inc[v] = m
                                for v, m in zip(vY, cy):
                                    if m:
                                        inc[v] = m
                                wgt = 1.0
                                for m in inc.values():
                                    wgt /= math.factorial(m)
                                for E in edge_multisets(slots, 1):
                                    links = [list(inc.keys())] + [list(k) for k in E]
                                    if not connected(nv, links):
                                        continue
                                    es = {}
                                    for (v, w), m in E.items():
                                        key = tuple(sorted((vsite[v], vsite[w])))
                                        es[key] = es.get(key, 0) + m
                                    units = {"D21": 2, "K31": 3, "K22": 2}[name]
                                    # random-signed hyperedge: cross-paired, identifies both of its sites
                                    es_r = dict(es); es_r[(min(X, Y), max(X, Y)) + ("h",)] = 1
                                    expo_r = exponent_h(legrows, P, vsite, ns, es, (X, Y), False) - units
                                    expo_c = exponent_h(legrows, P, vsite, ns, es, (X, Y), True) - units
                                    expo = max(expo_r, expo_c)
                                    yield P, vsite, ns, E, (name, X, Y, inc), wgt, expo


def _contract(site_ops, pair_ops, ns, outdiag):
    """Contract one structure (sites 0..ns-1 summed, rows i, j free) with BLAS matrix products only.
    site_ops[s] = {"i": n x n matrix or None, "j": ..., "v": vector}; pair_ops = list of (matrix, s, t).
    Leaf elimination: a site with at most one row letter and one neighbour is folded into the neighbour by one
    GEMM; a row-free site with two neighbours becomes a new edge; the last site gives A B^T (or a row sum).
    Returns None when no such order exists (caller falls back to einsum)."""
    S = {s: dict(site_ops[s]) for s in range(ns)}
    E = {}
    for M, s, t in pair_ops:
        key = (s, t) if s < t else (t, s)
        Mo = M if s < t else M.T
        E[key] = Mo if key not in E else E[key] * Mo
    alive = set(range(ns))
    def nbrs(s):
        return [t for t in alive if t != s and ((s, t) in E or (t, s) in E)]
    def edge(s, t):
        return E[(s, t)] if (s, t) in E else E[(t, s)].T     # oriented s -> t
    def mul_row(t, r, X):
        S[t][r] = X if S[t].get(r) is None else S[t][r] * X
    while len(alive) > 1:
        done = False
        for s in sorted(alive):
            nb = nbrs(s)
            rows = [r for r in ("i", "j") if S[s].get(r) is not None]
            v = S[s]["v"]
            if len(nb) == 1 and len(rows) <= 1:
                t = nb[0]; Est = edge(s, t)
                if rows:
                    r = rows[0]
                    X = (S[s][r] * v[None, :]) @ Est        # (r, t)
                    mul_row(t, r, X)
                else:
                    S[t]["v"] = S[t]["v"] * (v @ Est)
                E.pop((s, t) if (s, t) in E else (t, s)); alive.discard(s); done = True; break
            if len(nb) == 2 and not rows:
                t, u = nb
                X = (edge(t, s) * v[None, :]) @ edge(s, u)     # (t, u)
                E.pop((s, t) if (s, t) in E else (t, s)); E.pop((s, u) if (s, u) in E else (u, s))
                key = (t, u) if t < u else (u, t)
                Xo = X if t < u else X.T
                E[key] = Xo if key not in E else E[key] * Xo
                alive.discard(s); done = True; break
        if not done:
            return None
    s = alive.pop()
    A, B, v = S[s].get("i"), S[s].get("j"), S[s]["v"]
    if outdiag:
        return (A * v[None, :]).sum(1)
    if A is None or B is None:
        return None
    return (A * v[None, :]) @ B.T


def evaluate(legrows, st, F, mats, outdiag=False, pairs=None):
    """Exact sum over DISTINCT site indices for distinct site blocks (Moebius inversion over coincidences), each
    term contracted with BLAS products (einsum only as a fallback)."""
    P, vsite, ns, E, hed, wgt, expo = st
    deg = [0] * len(P)
    for (v, w), m in E.items():
        deg[v] += m; deg[w] += m
    if hed is not None:
        for v, m in hed[3].items():
            deg[v] += m
    letters = "abcd"
    site_raw = []
    for si in range(ns):
        vec = np.ones(n)
        pi = qj = 0
        for v in range(len(P)):
            if vsite[v] != si:
                continue
            vec = vec * F[(len(P[v]), deg[v])]
            for lg in P[v]:
                if legrows[lg] == "i":
                    pi += 1
                else:
                    qj += 1
        site_raw.append((pi, qj, vec))
    pair_ops = [(mats["Co"] ** m, vsite[v], vsite[w]) for (v, w), m in E.items()]
    if hed is not None:
        pair_ops.append((mats[hed[0]], hed[1], hed[2]))
    adjacent = {(min(x, y), max(x, y)) for _, x, y in pair_ops}
    total = 0.0
    for pi_ in set_partitions(range(ns)):
        if any((min(x, y), max(x, y)) in adjacent for blk in pi_ for x in blk for y in blk if x != y):
            continue
        rep = {}
        for bi, blk in enumerate(pi_):
            for x in blk:
                rep[x] = bi
        mob = 1
        for blk in pi_:
            mob *= (-1) ** (len(blk) - 1) * math.factorial(len(blk) - 1)
        nsm = len(pi_)
        sops = [{"i": None, "j": None, "v": np.ones(n)} for _ in range(nsm)]
        pcount = [[0, 0] for _ in range(nsm)]
        for si in range(ns):
            pi, qj, vec = site_raw[si]
            b = rep[si]
            sops[b]["v"] = sops[b]["v"] * vec
            pcount[b][0] += pi; pcount[b][1] += qj
        for b in range(nsm):
            if pcount[b][0]:
                sops[b]["i"] = mats["W"] ** pcount[b][0]
            if pcount[b][1]:
                sops[b]["j"] = mats["W"] ** pcount[b][1]
        if outdiag:
            for b in range(nsm):
                if sops[b]["j"] is not None:
                    sops[b]["i"] = sops[b]["j"] if sops[b]["i"] is None else sops[b]["i"] * sops[b]["j"]
                    sops[b]["j"] = None
        pops = [(M, rep[x], rep[y]) for M, x, y in pair_ops]
        if pairs is not None:
            for b in range(nsm):
                sops[b]["pi"], sops[b]["qj"] = pcount[b]
            val = _contract_pairs(sops, pops, nsm, pairs[0], pairs[1])
            if val is None:
                SKIPPED[0] += 1
                continue
            total = total + mob * val
            continue
        val = _contract(sops, pops, nsm, outdiag)
        if val is None and not EXACT_FALLBACK:
            SKIPPED[0] += 1
            continue
        if val is None:
            ops, subs = [], []
            for b in range(nsm):
                for r in ("i", "j"):
                    if sops[b][r] is not None:
                        ops.append(sops[b][r] * (sops[b]["v"][None, :] if r == "i" or sops[b]["i"] is None else 1.0))
                        subs.append(r + letters[b])
                if sops[b]["i"] is None and sops[b]["j"] is None:
                    ops.append(sops[b]["v"]); subs.append(letters[b])
            for M, x, y in pops:
                ops.append(M); subs.append(letters[x] + letters[y])
            expr = ",".join(subs) + "->" + ("i" if outdiag else "ij")
            FALLBACK[0] += 1
            val = np.einsum(expr, *ops, optimize="optimal")
        total = total + mob * val
    return wgt * total


FALLBACK = [0]


def _contract_pairs(site_ops, pair_ops, ns, I, J):
    """Pair-sampled contraction: rows i, j replaced by the sampled pairs (I[p], J[p]); every site operand is P x n.
    Tree elimination by GEMMs (P x n) @ (n x n). Returns None for cyclic site graphs."""
    X = {}
    for s in range(ns):
        op = site_ops[s]["v"][None, :] * np.ones((len(I), 1))
        if site_ops[s]["pi"]:
            op = op * mats_W[I] ** site_ops[s]["pi"]
        if site_ops[s]["qj"]:
            op = op * mats_W[J] ** site_ops[s]["qj"]
        X[s] = op
    E = {}
    for M, s, t in pair_ops:
        key = (s, t) if s < t else (t, s)
        Mo = M if s < t else M.T
        E[key] = Mo if key not in E else E[key] * Mo
    alive = set(range(ns))
    def nbrs(s):
        return [t for t in alive if t != s and ((s, t) in E or (t, s) in E)]
    while len(alive) > 1:
        leaf = next((s for s in sorted(alive) if len(nbrs(s)) == 1), None)
        if leaf is None:
            return None
        t = nbrs(leaf)[0]
        Est = E[(leaf, t)] if (leaf, t) in E else E[(t, leaf)].T
        X[t] = X[t] * (X[leaf] @ Est)
        E.pop((leaf, t) if (leaf, t) in E else (t, leaf)); alive.discard(leaf)
    return X[alive.pop()].sum(1)


mats_W = None
SKIPPED = [0]


def r2(Y, Tk, wmat, mask):
    T0, T1, Tf = Tk
    ip = lambda A, B: float(np.sum((A * B * wmat * wmat)[mask]))
    ET = ip(T0, T1)
    return 1.0 - ip(T0 - Y, T1 - Y) / ET, ip(Y, Tf) / ip(Y, Y)


lay = {h: {k: f64(T[h][k][L]) for k in ("mu", "var", "k3", "k4", "cov", "D21", "K22", "K31")} for h in ("full",)}["full"]
nxt = {h: {k: f64(T[h][k][L + 1]) for k in ("mu", "var", "k4", "K22", "K31", "cov")} for h in ("full", "h0", "h1")}
Co = 0.5 * (lay["cov"] + lay["cov"].T); np.fill_diagonal(Co, 0.0)
mats = {"W": W, "Co": Co}
for k in ("D21", "K31"):
    X = lay[k].copy(); np.fill_diagonal(X, 0.0); mats[k] = X
X = 0.5 * (lay["K22"] + lay["K22"].T); np.fill_diagonal(X, 0.0); mats["K22"] = X
Ft = site_factors(lay["mu"], lay["var"], lay["k3"], lay["k4"], True)
Fg = site_factors(lay["mu"], lay["var"], lay["k3"], lay["k4"], False)  # used only when GSF=1

mu1, v1 = nxt["full"]["mu"], nxt["full"]["var"]
s1 = np.sqrt(v1); a1 = mu1 / s1
act = a1 > -2.5
mask = off & act[:, None] & act[None, :]
ph1 = np.exp(-0.5 * a1 * a1) / SQ2PI
J1, J2, J3 = ndtr(a1), ph1 / s1, a1 * ph1 / s1**2       # J3 = E relu''' = -p'(0) = alpha phi / s^2
print(f"=== network {net}, transition {L} -> {L + 1}: {act.sum()} active neurons", flush=True)

mats_W = W
rng = np.random.default_rng(1000 * net + L)
ia, ja = np.nonzero(mask)
pick = rng.choice(len(ia), size=min(NPAIRS, len(ia)), replace=False)
PI, PJ = ia[pick], ja[pick]


def r2v(y, t0, t1, tf, w):
    ip = lambda A, B: float(np.sum(A * B * w * w))
    return 1.0 - ip(t0 - y, t1 - y) / ip(t0, t1), ip(y, tf) / ip(y, y)


for slice_name, legrows, tkey in (("(3,1)", "iiij", "K31"), ("(2,2)", "iijj", "K22"), ("diag", "iiii", "k4")):
    if {"(3,1)": "31", "(2,2)": "22", "diag": "diag"}[slice_name] not in SLICES:
        continue
    outdiag = slice_name == "diag"
    pairmode = slice_name == "(2,2)"
    if outdiag:
        sel = lambda X: X[act]
        w = np.ones(int(act.sum()))
        tv = [nxt[h]["k4"][act] for h in ("h0", "h1", "full")]
    else:
        if pairmode:
            sel = lambda X: X if X.ndim == 1 else X[PI, PJ]
            w = (J2[PI] * J2[PJ])
        else:
            sel = lambda X: X[mask]
            w = (J3[:, None] * J1[None, :])[mask]
        tv = []
        for h in ("h0", "h1", "full"):
            X = nxt[h][tkey]
            if tkey == "K22":
                X = 0.5 * (X + X.T)
            tv.append(X[PI, PJ] if pairmode else X[mask])
    groups = {}
    nst = 0
    SKIPPED[0] = 0
    for st in structures(legrows, hyper=True):
        expo = st[6]
        if expo < EXPO_MIN:
            continue
        nst += 1
        Y = evaluate(legrows, st, Ft, mats, outdiag, pairs=(PI, PJ) if pairmode else None)
        y = sel(Y) if not np.isscalar(Y) else np.zeros_like(w)
        P, vsite, ns, E, hed, _, _ = st
        e = sum(E.values())
        key = (("hyper:" + hed[0] + f" e={e} expo={expo}") if hed is not None else f"C e={e} expo={expo}")
        g = groups.setdefault(key, [0.0, 0])
        g[0] = g[0] + y; g[1] += 1
    ytot = sum(g[0] for g in groups.values())
    yc = sum(g[0] for k, g in groups.items() if not k.startswith("hyper"))
    R, sc = r2v(ytot, *tv, w)
    Rc, scc = r2v(yc, *tv, w)
    print(f"  {slice_name}: {nst} structures ({SKIPPED[0]} cyclic coincidence terms skipped) | all: R2 {100 * R:6.1f}% "
          f"scale {sc:+.3f} | C-only: R2 {100 * Rc:6.1f}% scale {scc:+.3f}"
          + (f" | on {len(PI)} sampled pairs" if pairmode else ""), flush=True)
    ipw = lambda A, B: float(np.sum(A * B * w * w))
    EY = ipw(ytot, ytot)
    for k in sorted(groups):
        g = groups[k]
        Rk, sk = r2v(g[0], *tv, w)
        print(f"      {k:24s} n={g[1]:3d}  share of prediction energy {100 * ipw(g[0], ytot) / EY:6.1f}%  "
              f"alone R2 {100 * Rk:6.1f}%  alone scale {sk:+.3f}", flush=True)
    if not outdiag:
        dilm = (v1[:, None] * nxt["full"]["cov"]) if tkey == "K31" else np.outer(v1, v1)
        dil = dilm[PI, PJ] if pairmode else dilm[mask]
        Rd, _ = r2v(dil * (ipw(dil, tv[2]) / ipw(dil, dil)), *tv, w)
        Rds, _ = r2v(ytot + dil * (ipw(dil, tv[2] - ytot) / ipw(dil, dil)), *tv, w)
        print(f"      dilation shape alone (best-scaled on the target): R2 {100 * Rd:6.1f}% | "
              f"all + best-scaled dilation shape: R2 {100 * Rds:6.1f}%", flush=True)
