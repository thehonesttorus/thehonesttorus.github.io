# The two-channel ledger in common / disagreement coordinates: g3 (common gain) and d = g4 - g3.  Exactly,
#   g3_{l+1} = g3_l + a d_l + s3 + e3,   d_{l+1} = (1 - a - b) d_l + (s4 - s3) + (e4 - e3),
# with a = r34, b = r43 of the row-normalised transport, s = one-loop generation, e = what the retained two-channel
# model misses (one-step residual against the measured gains of the next layer, no drift).
import numpy as np, sys
net = int(sys.argv[1]) if len(sys.argv) > 1 else 0
Z = np.load(f"gac3_off{net}.npz"); rows, gm3, gm4 = Z["rows"], Z["gm3"], Z["gm4"]
print(f"net {net}: common / disagreement ledger (one step from the measured gains of layer l)")
print("  l+1 | g3 meas  d meas | retention 1-a-b | source s4-s3 | transport of d (raw M) | residual e3  e4  e4-e3 | residual as % of g3, d")
for r in rows:
    l = int(r[0]); p3, p4 = r[3], r[8]; r33, r34, r43, r44 = r[11], r[12], r[13], r[14]
    g3, g4 = gm3[l], gm4[l]; d = g4 - g3
    pred3 = r33 * g3 + r34 * g4 + p3; pred4 = r43 * g3 + r44 * g4 + p4
    e3 = gm3[l + 1] - pred3; e4 = gm4[l + 1] - pred4
    a = r34 / (r33 + r34); b = r43 / (r43 + r44)
    dn = gm4[l + 1] - gm3[l + 1]
    print(f"  {l+1:3d} | {gm3[l+1]:+.5f} {dn:+.5f} | {1-a-b:.3f} | {p4-p3:+.5f} | {pred4-pred3:+.5f} | {e3:+.5f} {e4:+.5f} {e4-e3:+.5f} | {100*e3/gm3[l+1]:+.1f}% {100*(e4-e3)/abs(dn) if abs(dn)>1e-4 else float('nan'):+.0f}%")
