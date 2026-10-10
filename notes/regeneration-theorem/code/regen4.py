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
EXPO_MIN = int(os.environ.get("EXPO_MIN", "1"))   # Lemma 2 cut (1 = leading only; 0 adds the n^(-1/2) structures)   # also evaluate with Gaussian site factors (consistency check)
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


def evaluate(legrows, st, F, mats, outdiag=False):
    """Exact sum over DISTINCT site indices for distinct site blocks, by Moebius inversion over site coincidences."""
    P, vsite, ns, E, hed, wgt, expo = st
    deg = [0] * len(P)
    for (v, w), m in E.items():
        deg[v] += m; deg[w] += m
    if hed is not None:
        for v, m in hed[3].items():
            deg[v] += m
    letters = "abcd"
    site_ops = []          # per site: list of (operand, row-letter or None)
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
        ops = []
        if pi:
            ops.append((mats["W"] ** pi * vec[None, :], "i"))
            if qj:
                ops.append((mats["W"] ** qj, "j"))
        elif qj:
            ops.append((mats["W"] ** qj * vec[None, :], "j"))
        site_ops.append(ops)
    pair_ops = [(mats["Co"] ** m, vsite[v], vsite[w]) for (v, w), m in E.items()]
    if hed is not None:
        pair_ops.append((mats[hed[0]], hed[1], hed[2]))
    adjacent = {(min(x, y), max(x, y)) for _, x, y in pair_ops}
    out = "i" if outdiag else "ij"
    total = 0.0
    for pi_ in set_partitions(range(ns)):
        if any((min(x, y), max(x, y)) in adjacent for blk in pi_ for x in blk for y in blk if x != y):
            continue
        rep = {}
        for bi, blk in enumerate(pi_):
            for x in blk:
                rep[x] = letters[bi]
        mob = 1
        for blk in pi_:
            mob *= (-1) ** (len(blk) - 1) * math.factorial(len(blk) - 1)
        ops, subs = [], []
        for si in range(ns):
            for op, row in site_ops[si]:
                ops.append(op); subs.append(row + rep[si])
        for op, x, y in pair_ops:
            ops.append(op); subs.append(rep[x] + rep[y])
        expr = ",".join(subs) + "->" + out
        path = np.einsum_path(expr, *ops, optimize="greedy")
        big = float(re.search(r"Largest intermediate:\s*([0-9.eE+-]+)", path[1]).group(1))
        if big > 4 * n * n:
            path = np.einsum_path(expr, *ops, optimize="optimal")
        total = total + mob * np.einsum(expr, *ops, optimize=path[0])
    return wgt * total


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

for slice_name, legrows, wmat, tkey in (("(3,1)", "iiij", J3[:, None] * J1[None, :], "K31"),
                                        ("(2,2)", "iijj", J2[:, None] * J2[None, :], "K22"),
                                        ("diag", "iiii", None, "k4")):
    if {"(3,1)": "31", "(2,2)": "22", "diag": "diag"}[slice_name] not in SLICES:
        continue
    outdiag = slice_name == "diag"
    if outdiag:
        Tk = tuple(nxt[h]["k4"] for h in ("h0", "h1", "full"))
        msk = act; wm = np.ones(n)
    else:
        Tk = []
        for h in ("h0", "h1", "full"):
            X = nxt[h][tkey].copy()
            if tkey == "K22":
                X = 0.5 * (X + X.T)
            np.fill_diagonal(X, 0.0); Tk.append(X)
        Tk = tuple(Tk); msk = mask; wm = wmat
    groups = {}
    nst = 0
    for st in structures(legrows, hyper=True):
        expo = st[6]
        if expo < EXPO_MIN:
            continue
        nst += 1
        Yt = evaluate(legrows, st, Ft, mats, outdiag)
        Yg = evaluate(legrows, st, Fg, mats, outdiag) if GSF else 0.0 * Yt
        if not outdiag:
            np.fill_diagonal(Yt, 0.0); np.fill_diagonal(Yg, 0.0)
        P, vsite, ns, E, hed, _, _ = st
        e = sum(E.values())
        key = (("hyper:" + hed[0] + f" e={e} expo={expo}") if hed is not None else f"C e={e} expo={expo}")
        g = groups.setdefault(key, [0.0, 0.0, 0])
        g[0] = g[0] + Yt; g[1] = g[1] + Yg; g[2] += 1
    Ytot = sum(g[0] for g in groups.values())
    Ygtot = sum(g[1] for g in groups.values())
    Yc = sum(g[0] for k, g in groups.items() if not k.startswith("hyper"))
    R, sc = r2(Ytot, Tk, wm, msk)
    Rg, scg = r2(Ygtot, Tk, wm, msk) if GSF else (float("nan"), float("nan"))
    Rc, scc = r2(Yc, Tk, wm, msk)
    print(f"  {slice_name}: {nst} structures | all: R2 {100 * R:6.1f}% scale {sc:+.3f} | Gaussian site factors: "
          f"R2 {100 * Rg:6.1f}% scale {scg:+.3f} | C-only: R2 {100 * Rc:6.1f}% scale {scc:+.3f}", flush=True)
    ipw = lambda A, B: float(np.sum((A * B * wm * wm)[msk]))
    EY = ipw(Ytot, Ytot)
    for k in sorted(groups):
        g = groups[k]
        Rk, sk = r2(g[0], Tk, wm, msk)
        print(f"      {k:22s} n={g[2]:3d}  share of prediction energy {100 * ipw(g[0], Ytot) / EY:6.1f}%  "
              f"alone R2 {100 * Rk:6.1f}%", flush=True)
    if not outdiag:
        # production-like shapes, best-scaled (upper bound for a one-amplitude model)
        dil = (v1[:, None] * nxt["full"]["cov"]) if tkey == "K31" else np.outer(v1, v1)
        dil = dil.copy(); np.fill_diagonal(dil, 0.0)
        Rd, sd = r2(dil * (ipw(dil, Tk[2]) / ipw(dil, dil)), Tk, wm, msk)
        print(f"      dilation shape (best-scaled on the target): R2 {100 * Rd:6.1f}%", flush=True)
        Rds, _ = r2(Ytot + dil * (ipw(dil, Tk[2] - Ytot) / ipw(dil, dil)), Tk, wm, msk)
        print(f"      all + best-scaled dilation shape: R2 {100 * Rds:6.1f}%", flush=True)
