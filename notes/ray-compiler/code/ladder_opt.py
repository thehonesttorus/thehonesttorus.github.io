# Rank ladder from the readout law and the error-share law (ray-compiler note section 7.4), calibrated on two measurements.
#   python ladder_opt.py TAILS_FILE [TAILS_FILE ...]      (outputs of legspec.py NET MC --tails, averaged)
# Model: Delta MSE / MSE_dense = sum_a kappa_a eps(a, r_a), kappa_a = kappa_5 q^(a - 5), eps(a, r) the leg tails.
# Calibration (32 networks, paired against production = dense to age 4, 320 for ages 5-8, 192 nested for ages >= 9):
#   P1  all dense:            MSE_prod / MSE_dense - 1 = X1          (measured)
#   P2  one tier (all at 320): MSE_prod / MSE_onetier - 1 = X2       (measured)
# Cost model (units of 2 n^3 FLOPs, cold call): a confined member-layer costs G r (forming + factor-space D21 hub),
# a dense young member-layer D; each layer's join costs J r1 (range finder, factors, post-W projection); the nested
# move costs M per layer. G, D, J, M from the profile of note XXIX and the two measured FLOP differences.
import sys, re, math, numpy as np
ranks = (64, 96, 128, 160, 192, 224, 256, 288, 320, 352, 384, 448, 512)
tails = {}
files = [x for i, x in enumerate(sys.argv[1:], 1) if not x.startswith('--') and not sys.argv[i - 1].startswith('--')]
for f in files:
    for l in open(f):
        m = re.match(r"\s+a=\s*(\d+) \|.*?\| (.*?) \| free", l)
        if m:
            tails.setdefault(int(m.group(1)), []).append([float(x) for x in m.group(2).split()])
eps_tab = {a: np.mean(v, axis=0) for a, v in tails.items()}
amax = max(eps_tab)


def eps(a, r):
    """log-log interpolation of the tail in rank; ages beyond the measured range use the oldest measured age."""
    e = eps_tab[min(a, amax)]
    lr, le = np.log(ranks), np.log(np.maximum(e, 1e-12))
    return float(np.exp(np.interp(math.log(r), lr, le)))


L = 16
N = {a: (10 if a == 5 else 16 - a) for a in range(5, 16)}       # member-layers of age a in the old tiers (no join at 15)
N4 = 11                                                          # member-layers of age 4 if it were confined (layers 4-14)
X1, X2 = float(sys.argv[sys.argv.index("--x1") + 1]) if "--x1" in sys.argv else 0.0522, \
    float(sys.argv[sys.argv.index("--x2") + 1]) if "--x2" in sys.argv else 0.0045


def ladder_err(rk, k5, q):
    return sum(k5 * q ** (a - 5) * eps(a, rk(a)) for a in range(5, 16))


prod = lambda a: 320 if a <= 8 else 192
one = lambda a: 320
# solve for (kappa_5, q): X2 = err(prod) - err(one) relative to MSE_dense, X1 = err(prod)
best = None
for q in np.linspace(0.2, 1.0, 801):
    d = ladder_err(prod, 1.0, q) - ladder_err(one, 1.0, q)
    k5 = X1 / ladder_err(prod, 1.0, q)
    r = k5 * d
    if best is None or abs(r - X2 * (1 + X1)) < abs(best[2] - X2 * (1 + X1)):
        best = (q, k5, r)
q, k5, _ = best
print(f"calibrated: kappa_5 = {k5:.2f}, q = {q:.3f} (kappa halves every {math.log(0.5) / math.log(q):.2f} ages)")
print("  per-age share of the production confinement error (% of MSE_dense): "
      + " ".join(f"a{a}:{100 * k5 * q ** (a - 5) * eps(a, prod(a)):.2f}" for a in range(5, 16)))
G, J, D, Mv = 2.32e-3, 8.8e-3, 2.31, 1.5   # G from P2's FLOP difference, D from P1's, J from the profile
Cprod = 241.0                          # units incl. the residual-time term at production


def cost(rk, r1, nest=True):
    return G * sum(N[a] * rk(a) for a in N) + J * r1 * 10 + (Mv if nest else 0.0)


c_prod = cost(prod, 320)
for name, rk, r1, nest in (("production", prod, 320, True), ("one tier 320", one, 320, False)):
    e = ladder_err(rk, k5, q); c = cost(rk, r1, nest)
    print(f"  {name:28s} err {100 * e:5.2f}%  old-tier cost {c:6.1f} u  score vs production "
          f"{100 * ((1 + e) / (1 + ladder_err(prod, k5, q)) * (Cprod + c - c_prod) / Cprod - 1):+.2f}%")
# optimal ladders: 2 tiers (gate g2, ranks r1 >= r2) and 3 tiers (gates g2 < g3, ranks r1 >= r2 >= r3)
grid = range(160, 449, 16)
best2 = min(((ladder_err(lambda a: r1 if a < g2 else r2, k5, q), cost(lambda a: r1 if a < g2 else r2, r1), g2, r1, r2)
             for g2 in range(6, 13) for r1 in grid for r2 in range(64, r1 + 1, 16)),
            key=lambda t: (1 + t[0]) * (Cprod + t[1] - c_prod))
e0 = ladder_err(prod, k5, q)
f = lambda e, c: 100 * ((1 + e) / (1 + e0) * (Cprod + c - c_prod) / Cprod - 1)
print(f"  best 2-tier: ages 5-{best2[2] - 1} at {best2[3]}, ages >= {best2[2]} at {best2[4]}: err {100 * best2[0]:.2f}%, "
      f"score vs production {f(best2[0], best2[1]):+.2f}%")
best3 = None
for g2 in range(6, 11):
    for g3 in range(g2 + 1, 14):
        for r1 in range(256, 449, 32):
            for r2 in range(128, r1 + 1, 32):
                for r3 in range(64, r2 + 1, 32):
                    rk = lambda a, g2=g2, g3=g3, r1=r1, r2=r2, r3=r3: r1 if a < g2 else (r2 if a < g3 else r3)
                    e, c = ladder_err(rk, k5, q), cost(rk, r1) + Mv
                    s = (1 + e) * (Cprod + c - c_prod)
                    if best3 is None or s < best3[0]:
                        best3 = (s, e, c, g2, g3, r1, r2, r3)
_, e, c, g2, g3, r1, r2, r3 = best3
print(f"  best 3-tier: ages 5-{g2 - 1} at {r1}, {g2}-{g3 - 1} at {r2}, >= {g3} at {r3}: err {100 * e:.2f}%, "
      f"score vs production {f(e, c):+.2f}% (one extra nested move per layer charged)")
# the dense boundary: confine age 4 too (AGE_OLD = 3) in the same shared basis, one more join (layer 4)
for r1 in (320, 352, 384, 416, 448):
    rk = lambda a, r1=r1: r1 if a <= 8 else 192
    e = ladder_err(rk, k5, q) + k5 / q * eps(4, r1)
    c = cost(rk, r1) + N4 * (G * r1 - D) + J * r1
    print(f"  age 4 confined too, tier 1 at {r1}, nested 192: err {100 * e:.2f}%, score vs production {f(e, c):+.2f}%")
