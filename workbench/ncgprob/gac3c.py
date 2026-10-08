# The two-gain recursion with the fourth-cumulant channel completed by its dropped classes:
#   generation phi4_full(l) = one-loop kappa_4 of z_{l+1} by all five classes (k4classes_gauss.py; the pair ledger kept two),
#   transport r44_full(l) = r44_pair(l) + [6 sigma_diag^2 sigma_off^2 + 3 sigma_off^4] / (3 sigma^4) per unit gain (the mixture's
#   (2+1+1) and (1+1+1+1) terms, closed form), everything else as in gac3.py.  Variants: A = note XIX (pair M, pair phi);
#   B = pair M + mixture dropped-class transport, full phi4;  C = rows normalised, full phi4;  D = rows normalised, pair phi4.
import numpy as np, sys, os
net = int(sys.argv[1]) if len(sys.argv) > 1 else 0
Z = np.load(f"gac3_off{net}.npz"); rows, gm3, gm4 = Z["rows"], Z["gm3"], Z["gm4"]
G = {}
for l in range(1, 15):
    f = f"k4gauss_off{net}_l{l}.npz"
    if os.path.exists(f): G[l] = dict(np.load(f))
print(f"net {net}: closed two-gain recursions; one-loop by class available at source layers {sorted(G)}")
print("  l | measured g3 g4 ratio | A: XIX pair      | B: +dropped (mix transport, full phi4) | C: rows=1, full phi4 | D: rows=1, pair phi4 | phi4 pair -> full | r44 pair -> +mix")
st = {k: [0.0, 0.0] for k in "ABCD"}
for r in rows:
    l = int(r[0]); p3, p4 = r[3], r[8]; r33, r34, r43, r44 = r[11], r[12], r[13], r[14]
    if l in G: p4f = float(G[l]["gt"]); dr = float((G[l]["m211"] + G[l]["m1111"]) / gm4[l]) if gm4[l] > 1e-4 else 0.0
    else: p4f = p4; dr = 0.0
    s3, s4 = r33 + r34, r43 + r44
    M = {"A": (r33, r34, r43, r44, p3, p4), "B": (r33, r34, r43, r44 + dr, p3, p4f), "C": (r33 / s3, r34 / s3, r43 / s4, r44 / s4, p3, p4f), "D": (r33 / s3, r34 / s3, r43 / s4, r44 / s4, p3, p4)}
    for k, (a, b, c, d, q3, q4) in M.items():
        g3, g4 = st[k]; st[k] = [a * g3 + b * g4 + q3, c * g3 + d * g4 + q4]
    f = lambda k: f"{st[k][0]:.4f} {st[k][1]:.4f} {st[k][1]/st[k][0]:.2f}"
    print(f"  {l+1:2d} | {gm3[l+1]:.4f} {gm4[l+1]:.4f} {gm4[l+1]/gm3[l+1]:.2f} | {f('A')} | {f('B')} | {f('C')} | {f('D')} | {p4:+.5f} -> {p4f:+.5f} | {r44:.3f} -> {r44+dr:.3f}")
