"""Flopscope 0.12.1 cost calculator for closure-based K=3 chains at the Phase 2 shape (n = 1024, L = 16).

1 unit = one dense 1024^3 matmul = 2^31 FLOPs; B = 2^41 FLOPs = 1024 units.

Every price is a flopscope formula (matmul m w (2k-1), aliased Gram (2k-1) m(m+1)/2, Strassen-Winograd
recursion of the V29 kernel with its combo/assemble element ops, QR 2mnk-2k^3/3 x2, elementwise 1/elt,
transcendental 16/elt, stats.norm in float64) and is checked against the metered numbers in
ops_measured.json (probe_ops.py), designs_measured.json (probe_designs.py) and v29_ledger.json
(audit_v29.py: 504aldo's namespace-tagged V29 run under the meter on a random He MLP) by --check.

Designs (per-layer convention of the V29 ledger: layer l bills the transport l-1 -> l and the
nonlinearity at l; layer 0 bills the input Gram; layer 15 is trimmed to the mean, i.e. only diagonals):
  i    the published V29 structure (young dense sources, shared-basis old tiers, thin legs), Strassen L5
  ii   memoryless slices + identified closure births, only D21/D3 (and C) transported
  iii  ii + a rank-r (2,1,1) fourth-cumulant family (r symmetric modes, transported + born per layer)
  iv   ii + k transported covariance-response modes carrying old-source third-cumulant content
  v    K=2 covariance chain baseline

Usage:
  python cost.py                         # design table at the default (cheapest shippable) forms
  python cost.py --layers                # + per-layer units for every design
  python cost.py --check                 # formulas vs metered numbers (needs the *.json next to this file)
  python cost.py --lev 3 --r 8 --k 8 --legs C_phi,C_w2,D_w2,DT_phi --cov sym3 --births 1
"""
import argparse
import json
import math
import os

N = 1024
L = 16
UNIT = float(2 ** 31)
B_UNITS = 1024.0
HERE = os.path.dirname(os.path.abspath(__file__))
STRASSEN_MIN = 32
FUSE_P = 343

# residual per flopscope call (ms), measured: V29 on this 4-core shared box 0.546 s / 13,121 calls
# = 0.042 ms; design skeletons 0.036-0.049 ms. V29 on its author's box (published) 0.25-0.32 s
# for the same 13,121 calls = 0.019-0.024 ms. Cap 400 ms per MLP; plan margin 2x -> 200 ms.
RESID_MS_PER_CALL_HERE = 0.042
RESID_MS_PER_CALL_GRADER = 0.022


# ----------------------------------------------------------------------------- primitive prices ---
def mm(m, k, w, b=1, rate=1):
    """flopscope matmul / einsum 'ij,jk->ik': m w (2k-1) per product."""
    return b * m * w * (2 * k - 1) * rate / UNIT


def gram_alias(rows, cols, rate=1):
    """einsum('ji,jk->ik', X, X) with the same object: (2 rows - 1) cols (cols+1)/2 (symmetric output)."""
    return (2 * rows - 1) * cols * (cols + 1) / 2 * rate / UNIT


def ew(numel, weight=1, rate=1):
    return numel * weight * rate / UNIT


def qr_reduced(m, k):
    return 2 * (2 * m * k * k - 2 * k ** 3 // 3) / UNIT


def _ok(m, kd, w, lev):
    d = 2 ** lev
    return lev > 0 and m % d == 0 and kd % d == 0 and w % d == 0 and min(m, kd, w) // d >= STRASSEN_MIN


def strassen_level(m, kd, w, lev):
    while lev > 0 and not _ok(m, kd, w, lev):
        lev -= 1
    return lev


def strassen(m, kd, w, lev, bx=1, by=1, P=1):
    """FLOPs (units) and flopscope calls of _Strassen.mm (504aldo V29) on X (bx,P,m,kd) @ Y (by,P,kd,w).
    Non-fused level: 7 left combos (5 add/sub + 2 copies, 1/elt) on bx P h q, 7 right combos on
    by P q v, recursion on (bx, by, 7P), 8 assemble adds on by P h v. Fused leaf (P >= 343 and the next
    level not allowed): 5 left ops, 5 right ops, 7 products, 8 accumulate ops."""
    if not _ok(m, kd, w, lev):
        return by * P * m * w * (2 * kd - 1) / UNIT, 1
    h, q, v = m // 2, kd // 2, w // 2
    if P >= FUSE_P and not _ok(h, q, v, lev - 1):
        f = 5 * bx * P * h * q + 5 * by * P * q * v + 7 * by * P * h * v * (2 * q - 1) + 8 * by * P * h * v
        return f / UNIT, 5 + 5 + 7 + 8
    sub, calls = strassen(h, q, v, lev - 1, bx, by, 7 * P)
    f = 7 * bx * P * h * q + 7 * by * P * q * v + 8 * by * P * h * v
    return f / UNIT + sub, calls + 7 + 7 + 8


def strassen_pool_bytes(m, kd, w, lev, bx=1, by=1, P=1, itemsize=4):
    """bytes of the pooled scratch the V29 kernel allocates for one family (combo buffers Xc, Yc and the
    product buffer M per level, which reuses Xc when the callee recurses on square blocks; fused leaf:
    one left, one right, one product block). Pools are keyed by geometry and shared by role across
    families of the same geometry (one _Strassen instance), so a design's peak is set by its largest
    family per geometry."""
    if not _ok(m, kd, w, lev):
        return 0
    h, q, v = m // 2, kd // 2, w // 2
    if P >= FUSE_P and not _ok(h, q, v, lev - 1):
        return itemsize * P * (bx * h * q + by * q * v + by * h * v)
    callee_recurses = _ok(h // 2, q // 2, v // 2, lev - 2)
    own = bx * 7 * P * h * q + by * 7 * P * q * v + (0 if (callee_recurses and q == v) else by * 7 * P * h * v)
    return itemsize * own + strassen_pool_bytes(h, q, v, lev - 1, bx, by, 7 * P, itemsize)


def strassen_hub(m, kd, w, lev, k, P=1):
    """_Strassen.hub: out (P,m,w) = sum_k X_k Y_k^T; dense fallback = batched GEMM + k-sum."""
    if not _ok(m, kd, w, lev):
        return (k * P * m * w * (2 * kd - 1) + (k - 1) * P * m * w) / UNIT, 2
    h, q, v = m // 2, kd // 2, w // 2
    if P >= FUSE_P and not _ok(h, q, v, lev - 1):
        f = 5 * k * P * h * q + 5 * k * P * v * q + 7 * (k * P * h * v * (2 * q - 1) + (k - 1) * P * h * v) + 8 * P * h * v
        return f / UNIT, 10 + 14 + 8
    sub, calls = strassen_hub(h, q, v, lev - 1, k, 7 * P)
    f = 7 * k * P * h * q + 7 * k * P * v * q + 8 * P * h * v
    return f / UNIT + sub, calls + 7 + 7 + 8


def family(b, lev, m=N, kd=N, w=N, shared_left=False):
    """b independent products (m,kd)@(kd,w) in one batched family (pooled operands, one op stream).
    Returns (units, calls) including the 2b operand loads (copyto, 1/elt) an estimator pays to stack
    the operands (shared_left: one left operand for the whole family, loaded once)."""
    lev = strassen_level(m, kd, w, lev)
    u, c = strassen(m, kd, w, lev, bx=1 if shared_left else b, by=b)
    loads = (1 if shared_left else b) * m * kd + b * kd * w
    return u + loads / UNIT, c + (1 if shared_left else b) + b


# level of the 3 half-size output blocks of a symmetric product: the kernel's own fallback
# (level(n/2, n, n/2, lev): L5 -> L4); probe_designs.py used lev - 1, which --check passes explicitly
BLEV_MODE = {"mode": "kernel"}


def blev(lev):
    if BLEV_MODE["mode"] == "minus1":
        return max(lev - 1, 0)
    return strassen_level(N // 2, N, N // 2, lev)


# ----------------------------------------------------------------------------- sandwich forms -----
def sandwich(form, lev):
    """W^T diag(p) S diag(p) W for symmetric S (n x n). Returns (units, calls, note).
    einsum_tagged: as_symmetric(S) + einsum('ji,jk,kl->il', X, S, X), X the same object: 1.5 u
      (HAZARD: raises SymmetryError 'max deviation = inf' on an indefinite symmetric S with O(1)
      entries, measured; fine on SPD inputs and small-scale ones)
    dense2: X.T @ (S @ X); strassen2: both products through Strassen
    sym3: T = S X (Strassen, can ride in another family as one slot) + the 3 independent (n/2, n) @ (n, n/2)
      output blocks one level shallower (V29's C_pre form)
    factor: S = U U^T available: Y = U^T X (Strassen) + aliased Gram Y^T Y (0.5 u); + Cholesky n^3/3 if not."""
    if form == "einsum_tagged":
        # stage 1 S X dense, stage 2 X^T T with symmetric output (orbit count), + the as_symmetric tag (7n^2 - 1)
        return mm(N, N, N) + gram_alias(N, N) + ew(7 * N * N - 1), 3, "SymmetryError hazard on indefinite O(1) S"
    if form == "dense2":
        return 2 * mm(N, N, N) + ew(N * N), 3, ""
    if form == "strassen2":
        u, c = strassen(N, N, N, strassen_level(N, N, N, lev))
        return 2 * u + ew(N * N), 2 * c + 1, ""
    if form == "sym3":
        u1, c1 = family(1, lev)
        u3, c3 = family(3, blev(lev), N // 2, N, N // 2)
        return u1 + u3 + ew(N * N) + ew(N * N), c1 + c3 + 6, ""
    if form == "factor":
        u, c = strassen(N, N, N, strassen_level(N, N, N, lev))
        return u + gram_alias(N, N) + ew(N * N), c + 2, "needs S = U U^T (Cholesky +0.167 u otherwise)"
    raise ValueError(form)


# ----------------------------------------------------------------------------- closure designs ----
# leg types of the identified first-order closure (notes/experiments/oracle_k3.py residual_basis):
# a star diagram with centre c and legs (E, x) costs, per distinct leg type, one a-leg product
# G[E,x] = E^T diag(x) W and one b-leg product E @ Lt (Lt = sum of vertex-weighted W * G over every
# term sharing that leg); every b-index-on-W contribution merges into ONE final L^T W.
LEGS = {
    "C_phi": "C_off with leaf weight Phi: leading Wick (centre w2), B2, B3 (D3 + two edges), B5 C-leg, B6 r=1 core",
    "C_w2": "C_off with leaf weight w2: B5 (K22 = lam C_off regenerated) K22-leg",
    "D_w2": "D21z^T with leaf weight w2: B1 (D21 hyperedge + edge j-k)",
    "DT_phi": "D21z with leaf weight Phi: B2 (D21 hyperedge + edge i-j)",
    "CC_x": "C_off o C_off leg: Gaussian rho^3 path terms (B4 without its n^4 triangle)",
}
DEFAULT_LEGS = ["C_phi", "C_w2", "D_w2", "DT_phi"]


def transition(legs, lev, cov="sym3", nmodes=0, births=1, last=False, first=False, slices_in_A=True,
               mode_basis=0, closure=True):
    """one transition s -> s+1 (s = l-1) of a closure chain, grouped into batched families the way an
    estimator would run it. Returns a list of (label, units, calls).

    A-family (one Strassen family): the a-leg products G[E,x] = E^T diag(x) W of every leg type, the
      covariance slot C_a W (cov = 'sym3'), the mode slots M^m (Phi W), and (slices_in_A) the two slice
      products S21 W and S21^T (W o W)
    B-family: the b-leg products E @ Lt, one per leg type (Lt = sum over the terms sharing the leg)
    F: the single final L^T W into which every b-index-on-W contribution merges
    3-block family: the 3 independent (n/2, n) @ (n, n/2) output blocks of each symmetric product
      (covariance + modes), one level shallower than the kernel's top level
    births: per mode, `births` weighted Grams G^T diag(d) G (split-sign aliased einsum, 0.503 u, dense)
    mode_basis q > 0: modes kept as U C_m U^T in a shared q-column basis instead of dense n x n:
      transport U (one (n,n)@(n,q)), project births (G_C U: one (n,n)@(n,q) + (n,q) Grams), refresh the basis
      (range finder: two (n,n)@(n,q) + QR(n,q)); r cores are q x q
    last (trimmed final layer): only diagonals are needed: a-legs + S21 W + C_a W, n^2 dot products."""
    legs = ["C_phi"] if first else legs
    nl = len(legs) if closure else 0
    out = []
    dense_modes = nmodes if mode_basis == 0 else 0
    if last:
        slots = nl + (1 if closure else 0) + 1 + dense_modes
        out.append((f"A-family ({slots} slots)", *family(slots, lev)))
        out.append(("n^2 diagonals", ew((4 * nl + 10 + 2 * dense_modes) * N * N), 4 * nl + 10 + 2 * dense_modes))
        if mode_basis:
            out.append(("mode basis diag", mm(N, N, mode_basis) + ew(nmodes * N * mode_basis * 2), 4))
        return out
    slots = nl + (1 if cov == "sym3" else 0) + dense_modes + (2 if (closure and slices_in_A) else 0)
    if slots:
        out.append((f"A-family ({slots} slots)", *family(slots, lev)))
    if closure and not slices_in_A:
        out.append(("slice family (2)", *family(2, lev)))
    if cov != "sym3":
        u, c, _ = sandwich(cov, lev)
        out.append(("cov sandwich " + cov, u, c))
    if closure:
        out.append((f"B-family ({nl})", *family(nl, lev)))
        out.append(("final L^T W", *family(1, lev)))
        out.append(("n^2 assembly", ew((6 * nl + 12) * N * N), 6 * nl + 12))
    nsym = (1 if cov == "sym3" else 0) + dense_modes
    if nsym:
        u3, c3 = family(3 * nsym, blev(lev), N // 2, N, N // 2)
        out.append((f"3-block family ({3 * nsym})", u3 + ew(2 * nsym * N * N), c3 + 2 * nsym))
    if nmodes and births:
        if mode_basis:
            q = mode_basis
            out.append(("mode births in basis", births * (mm(N, N, q) + nmodes * gram_alias(N, q)) + ew(nmodes * N * N), 4 + 3 * nmodes))
        else:
            out.append(("mode births (split-sign Gram)", nmodes * births * (gram_alias(N, N) + ew(6 * N * N)), nmodes * births * 9))
    if nmodes:
        if mode_basis:
            q = mode_basis
            out.append(("mode basis transport + refresh", 3 * mm(N, N, q) + qr_reduced(N, q) + nmodes * mm(q, q, q) * 2, 8 + 2 * nmodes))
            out.append(("mode D21 terms (factored)", nmodes * (mm(N, q, q) + ew(4 * N * q)) + mm(N, q, N), 3 * nmodes + 1))
        else:
            out.append(("mode vectors + D21 terms", ew(nmodes * 8 * N * N), 4 * nmodes))
    return out


def nonlin_layer(first=False):
    """the ReLU closure at one layer: Gaussian weights (norm.cdf/pdf in float64 on n-vectors, ~0.0001 u)
    and the term program on n x n (V29 ledger: nonlin 0.21 u / 70 calls per layer, layer 0 0.07 u)."""
    if first:
        return [("nonlin (layer 0)", 0.07, 160)]
    return [("nonlin + term program", 0.21, 70)]


def design_layers(design, lev=5, legs=None, cov="sym3", r=8, k=8, births=1, mode_basis=0):
    """per-layer list of (label, units, calls) for one design."""
    legs = DEFAULT_LEGS if legs is None else legs
    layers = []
    for li in range(L):
        if li == 0:
            layers.append([("input Gram W0^T W0 (aliased)", gram_alias(N, N), 1)] + nonlin_layer(first=True))
            continue
        nmodes = r if design == "iii" else (k if design == "iv" else 0)
        items = transition(legs, lev, cov, nmodes=nmodes, births=births, last=(li == L - 1), first=(li == 1),
                           mode_basis=mode_basis, closure=(design != "v"))
        layers.append(items + nonlin_layer())
    return layers


# ----------------------------------------------------------------------------- V29 replay ---------
def v29_layers(lev=5, age_old=4, r1=384, r2=224, age_old2=7):
    """structural replay of 504aldo's V29 op stream (estimator_v29.py; shapes from its metered op log).
    Dense n^3 families are priced with the Strassen formula at the kernel's own levels; the old tiers by
    their (n x r) shapes; the thin families by their (n x 16/18/32) shapes. Returns per-layer items."""
    layers = []
    for li in range(L):
        it = []
        k = li                                  # live sources at the transport into layer li
        last = li == L - 1
        if li == 0:
            it.append(("cpre Gram", gram_alias(N, N), 1))
            it.append(("birth", 0.07, 160))
            layers.append(it)
            continue
        # young family: W @ [A, P legs of the young sources, newborn A, C]; left operand W shared
        ny = min(k - 1, age_old - 1) if not last else age_old
        slots = 2 * ny + 2
        u, c = strassen(N, N, N, lev, bx=1, by=slots)
        it.append(("young_transport", u + ew(slots * N * N), c + slots))
        if not last:
            kh = 2 * min(k, age_old)
            u, c = strassen_hub(N, N, N, lev, kh)
            it.append(("hub", u + ew(2 * kh * N * N), c + kh))
            u3, c3 = strassen(N // 2, N, N // 2, strassen_level(N // 2, N, N // 2, lev), bx=3, by=3)
            it.append(("cpre (3 blocks)", u3 + ew(4 * N * N), c3 + 22))
        # old tiers: sources older than age_old joined to a shared rank-r1 basis (one join per layer
        # 5..14), those older than age_old2 nested in a rank-r2 sub-basis
        lj = li if not last else li - 1          # no join and no tier-2 move at the trimmed last layer
        n_old = max(0, lj - age_old)
        t2 = max(0, lj - age_old2)
        t1 = n_old - t2
        if not last and li > age_old:
            # range finder: 4 (n,n)@(n,r1) sketch products + QR; projections 2 (r1,n)@(n,n); Qc = W Qn;
            # Gram core 2 (r1,n)@(n,r1); from the second join on, the rotation of the tier-1 factors
            # (2 legs x (r1,r1)@(r1,n) per rotated source) and the small core updates
            join = 4 * mm(N, N, r1) + qr_reduced(N, r1) + mm(N, N, r1) + 2 * mm(r1, N, N) + 2 * mm(r1, N, r1)
            if li > age_old + 1:
                rot = min(li - age_old - 1, age_old2 - age_old)
                join += mm(r1, N, r1) + mm(N, r1, r1) + mm(r1, r1, r1) + mm(r1, N, r1) + 2 * mm(r1, r1, r1) + 2 * mm(r1, r1, N) * rot
            if li > age_old2 + 1:
                join += 0.049
            it.append(("join (range finder, projections, core, Qc, rotation)", join, 30))
            if t2 > 0:
                it.append(("tier-2 move", 2 * mm(r1, N, r1) + 2 * mm(r2, r1, N) + 3 * mm(r1, r1, r2) + 2 * qr_reduced(r1, r2) / 2
                           + 2 * mm(r2, r2, N) * (t2 - 1), 14))
        if n_old > 0:
            l1 = strassen_level(N, r1, N, lev)
            l2 = strassen_level(N, r2, N, lev)
            u1 = strassen(N, r1, N, l1)[0]
            u2 = strassen(N, r2, N, l2)[0]
            it.append(("old_legs (dense legs from factors)", 2 * t1 * u1 + 2 * t2 * u2, 30 + 10 * n_old))
            # (tier-1 legs (n,r1)@(r1,n) at Strassen L3, tier-2 (n,r2)@(r2,n) at L2)
            if not last:
                it.append(("shared (tier contractions + inner Qc^T)", 2 * t1 * u1 + 2 * t2 * u2 + mm(N, r1, N), 40 + 10 * n_old))
            else:
                it.append(("Qc transport", mm(N, N, r1), 1))
        # thin legs per live source: residual Z (n x 18) + feedback (n x 32), their contractions
        it.append(("thin_transport", k * mm(N, N, 50) if not last else (k + 1) * mm(N, N, 50) - mm(N, N, 50) + mm(N, N, 50), 4))
        # per-source thin families, constants read off the metered V29 ledger (L14 / 14 sources):
        # fb (rank-16 D21 feedback legs) 0.104, thin (rank-18 residual legs) 0.035, elem 0.028 u
        if not last:
            it.append(("thin + fb + elem (ledger per-source constants)", k * (0.104 + 0.0352 + 0.0279), 20 + 6 * k))
            it.append(("nonlin + birth + feed (ledger constants)", 0.40 + 0.0044 * k, 160))
        else:
            it.append(("fb + elem (D3 mode, ledger constants)", 0.535 + 0.352 + 0.029, 300))
        layers.append(it)
    return layers


# ----------------------------------------------------------------------------- reporting ----------
def design_memory_gb(design, lev=5, legs=None, cov="sym3", r=8, k=8, births=1, mode_basis=0):
    """peak memory estimate (GB) of one middle transition: the largest n x n Strassen family's pools
    (A-family: legs + cov slot + dense mode slots + 2 slices), the 3-block family's pools, the family
    operand/output stacks (3 n^2 per slot), and the persistent state (W stack, C, C_a, D21, S21, legs,
    modes). The real V29 peaks at 5.5 GB (its author's measurement); the grader cap is 8 GB."""
    legs = DEFAULT_LEGS if legs is None else legs
    nl = len(legs) if design != "v" else 0
    nm = (r if design == "iii" else (k if design == "iv" else 0)) if mode_basis == 0 else 0
    slots = nl + (1 if cov == "sym3" else 0) + nm + (2 if design != "v" else 0)
    l1 = strassen_level(N, N, N, lev)
    pool = strassen_pool_bytes(N, N, N, l1, slots, slots) if slots else 0
    nsym = (1 if cov == "sym3" else 0) + nm
    l3 = blev(lev)
    pool3 = strassen_pool_bytes(N // 2, N, N // 2, strassen_level(N // 2, N, N // 2, l3), 3 * nsym, 3 * nsym) if nsym else 0
    stacks = 4 * N * N * (3 * slots + 3 * nsym)
    state = 4 * N * N * (16 + 6 + 3 * nl + nm + (nm if births else 0)) + (4 * N * mode_basis * 3 if mode_basis else 0)
    return (pool + pool3 + stacks + state) / 1e9


def total(layers):
    per = [sum(x[1] for x in it) for it in layers]
    calls = [sum(x[2] for x in it) for it in layers]
    return per, calls


def table(configs):
    rows = []
    try:
        led = json.load(open(os.path.join(HERE, "v29_ledger.json")))
        v29_calls = led["call2"]["ops"]
    except (OSError, KeyError):
        v29_calls = 13121
    for name, layers in configs:
        per, calls = total(layers)
        T, C = sum(per), sum(calls)
        if name.startswith("i   V29"):
            C = v29_calls        # metered op count of the real V29 (the replay prices FLOPs, not calls)
        rows.append((name, per, T, T / B_UNITS, C, C * RESID_MS_PER_CALL_HERE, C * RESID_MS_PER_CALL_GRADER, MEM.get(name, float("nan"))))
    return rows


MEM = {}


def default_configs(lev=5, r=8, k=8, legs=None, cov="sym3", births=1, q=0):
    names = {"ii": f"ii  closure, {len(legs or DEFAULT_LEGS)} leg types, L{lev}", "iii": f"iii ii + (2,1,1) family r={r}, L{lev}",
             "iv": f"iv  ii + k={k} cov-response modes, L{lev}", "v": f"v   K=2 covariance, L{lev}",
             "iii-sb": f"iii-sb r={r} modes in a q={q} basis, L{lev}", "iv-sb": f"iv-sb k={k} modes in a q={q} basis, L{lev}"}
    MEM["i   V29 (replay, Strassen L5)"] = 5.5
    for key, nm in names.items():
        d = key.split("-")[0]
        MEM[nm] = design_memory_gb(d, lev, legs, cov, r, k, births, q if key.endswith("sb") else 0)
    cfg = [("i   V29 (replay, Strassen L5)", v29_layers(5)),
           (f"ii  closure, {len(legs or DEFAULT_LEGS)} leg types, L{lev}", design_layers("ii", lev, legs, cov)),
           (f"iii ii + (2,1,1) family r={r}, L{lev}", design_layers("iii", lev, legs, cov, r=r, births=births)),
           (f"iv  ii + k={k} cov-response modes, L{lev}", design_layers("iv", lev, legs, cov, k=k, births=births)),
           (f"v   K=2 covariance, L{lev}", design_layers("v", lev, legs, cov))]
    if q:
        cfg.insert(3, (f"iii-sb r={r} modes in a q={q} basis, L{lev}", design_layers("iii", lev, legs, cov, r=r, births=births, mode_basis=q)))
        cfg.insert(5, (f"iv-sb k={k} modes in a q={q} basis, L{lev}", design_layers("iv", lev, legs, cov, k=k, births=births, mode_basis=q)))
    return cfg


def check():
    ok = True
    ops = json.load(open(os.path.join(HERE, "ops_measured.json")))["results"]
    m = {r["name"].split(":")[0]: r for r in ops if "units" in r}
    print("== Strassen formula vs metered single product / batch / hub ==")
    for lev in range(6):
        f = strassen(N, N, N, lev)
        meas = [r for r in ops if r["name"].startswith(f"strassen L{lev}: one")][0]
        print(f"  L{lev}: formula {f[0]:.4f} u, {f[1]} calls | metered {meas['units']:.4f} u, {meas['calls']} calls")
        ok &= abs(f[0] - meas["units"]) < 2e-3
    for lev in (3, 5):
        f = strassen(N, N, N, lev, bx=4, by=4)
        meas = [r for r in ops if r["name"].startswith(f"strassen L{lev}: batch")][0]
        print(f"  L{lev} batch 4: formula {f[0]:.4f} | metered {meas['units']:.4f}")
        fh = strassen_hub(N, N, N, lev, 4)
        mh = [r for r in ops if r["name"].startswith(f"strassen hub L{lev}")][0]
        print(f"  L{lev} hub k=4: formula {fh[0]:.4f} | metered {mh['units']:.4f}")
    print("== sandwich forms ==")
    for form, lev, key in (("dense2", 0, "sandwich_plain"), ("einsum_tagged", 0, "sandwich_einsum3"),
                           ("strassen2", 5, "sandwich_strassen L5"), ("sym3", 5, "sandwich_strassen_sym3 L5"),
                           ("sym3", 4, "sandwich_strassen_sym3 L4"), ("factor", 5, "sandwich_factor_strassen L5")):
        rr = [r for r in ops if r["name"].startswith(key)][0]
        print(f"  {form:14s} L{lev}: formula {sandwich(form, lev)[0]:.4f} | metered {rr['units']:.4f} ({rr['calls']} calls)")
    try:
        des = json.load(open(os.path.join(HERE, "designs_measured.json")))["results"]
        print("== one middle transition: calculator vs metered skeleton ==")
        for r in des:
            if "design" not in r or "units" not in r:
                continue
            d, lev = r["design"], r["lev"]
            BLEV_MODE["mode"] = "minus1"     # the skeletons run the 3 output blocks at lev - 1
            legs = ["C_phi"] if d == "ii-wick" else DEFAULT_LEGS
            nm = r.get("r", 0) or r.get("k", 0)
            cv = "einsum_tagged" if r["cov"] == "einsum" else r["cov"]
            items = transition(legs, lev, cv, nmodes=nm, births=1 if d == "iv" else 0,
                               slices_in_A=False, closure=(d != "v"))
            BLEV_MODE["mode"] = "kernel"
            u = sum(x[1] for x in items)
            c = sum(x[2] for x in items)
            print(f"  {r['name'][:56]:56s} calc {u:7.3f} u {c:5d} calls | metered {r['units']:7.3f} u {r['calls']:5d} calls"
                  f" {r['residual_ms']:6.1f} ms")
    except FileNotFoundError:
        print("  (designs_measured.json not present)")
    led = json.load(open(os.path.join(HERE, "v29_ledger.json")))
    print(f"== V29 replay vs metered ledger (metered total {led['call2']['units']:.2f} u, {led['call2']['ops']} calls,"
          f" residual {led['call2']['resid']:.3f} s on this box) ==")
    per, calls = total(v29_layers())
    for li in range(L):
        lf = led["layer_family"][f"L{li:02d}"]
        lc = led["layer_family_calls"][f"L{li:02d}"]
        print(f"  L{li:02d}: replay {per[li]:6.2f} u | metered {sum(lf.values()):6.2f} u   calls replay {calls[li]:5d} | metered {sum(lc.values()):5d}")
    print(f"  total: replay {sum(per):.2f} u | metered {led['call2']['units']:.2f} u")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lev", type=int, default=5)
    ap.add_argument("--r", type=int, default=8)
    ap.add_argument("--k", type=int, default=8)
    ap.add_argument("--legs", default=",".join(DEFAULT_LEGS))
    ap.add_argument("--cov", default="sym3", help="sym3 | einsum (= einsum_tagged) | strassen2 | factor | dense2")
    ap.add_argument("--births", type=int, default=1)
    ap.add_argument("--q", type=int, default=0, help="also price iii/iv with modes in a shared q-column basis")
    ap.add_argument("--layers", action="store_true")
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    if a.check:
        check()
        return
    legs = [x for x in a.legs.split(",") if x]
    if a.cov == "einsum":          # the skeleton/REPORT name of the tagged 3-operand einsum covariance
        a.cov = "einsum_tagged"
    rows = table(default_configs(a.lev, a.r, a.k, legs, a.cov, a.births, a.q))
    print(f"{'design':44s} {'units':>8s} {'C/B':>7s} {'calls':>7s} {'resid here':>11s} {'resid grader':>13s} {'mem GB':>7s}")
    for name, per, T, cb, C, rh, rg, mem in rows:
        print(f"{name:44s} {T:8.1f} {cb:7.3f} {C:7d} {rh:9.0f}ms {rg:11.0f}ms {mem:7.1f}")
    if a.layers:
        print("\nper layer (units):")
        print(f"{'':44s} " + " ".join(f"L{li:02d}" for li in range(L)))
        for name, per, *_ in rows:
            print(f"{name:44s} " + " ".join(f"{x:4.1f}"[-4:] if x < 10 else f"{x:4.0f}" for x in per))


if __name__ == "__main__":
    main()
