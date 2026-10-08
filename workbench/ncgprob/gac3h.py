# Homogeneity-corrected two-gain recursion from the saved ledger: positive homogeneity (relu(G y) = G relu(y)) makes the
# scale mixture an exact invariant family, so the full transport of a unit-gain mixture is exactly 1 in both channels;
# the ledger's row sums (pair slices, fourth-order truncation) fall short.  Rescale each row to sum 1, keep the mixing
# ratio, and run (g3,g4)_{l+1} = M^_l (g3,g4)_l + (phi3,phi4)_l.  Also: the second eigenvalue of M^ (the relaxation rate
# of the non-mixture part of the fresh non-Gaussianity), and the stationary ratio it implies.
import numpy as np, sys
net = int(sys.argv[1]) if len(sys.argv) > 1 else 0
Z = np.load(f"gac3_off{net}.npz"); rows, gm3, gm4 = Z["rows"], Z["gm3"], Z["gm4"]
print(f"net {net}: homogeneity-corrected recursion (rows of M rescaled to sum 1)")
print("  l | row sums r3 r4 | mixing a34 a43 | lambda2 | measured g3 g4 ratio | derived g3 g4 ratio | measured-M derived ratio")
g3 = g4 = 0.0; h3 = h4 = 0.0
for r in rows:
    l = int(r[0]); p3, p4 = r[3], r[8]; r33, r34, r43, r44 = r[11], r[12], r[13], r[14]
    s3, s4 = r33 + r34, r43 + r44; a34 = r34 / s3; a43 = r43 / s4; lam2 = 1 - a34 - a43
    g3, g4 = (1 - a34) * g3 + a34 * g4 + p3, a43 * g3 + (1 - a43) * g4 + p4
    h3, h4 = r33 * h3 + r34 * h4 + p3, r43 * h3 + r44 * h4 + p4
    print(f"  {l+1:2d} | {s3:.3f} {s4:.3f} | {a34:.3f} {a43:.3f} | {lam2:.3f} | {gm3[l+1]:.5f} {gm4[l+1]:.5f} {gm4[l+1]/gm3[l+1]:.3f} | {g3:.5f} {g4:.5f} {g4/g3:.3f} | {h4/h3:.3f}")
