# Note XLIV section 9: pin targets from a recorded baseline run (V62_GPK=...:rec output lines "[gpkrec] layer L e_v ... e_31 ...").
#   python -I pin_from_rec.py REC.txt OUT.npz
import sys, re, numpy as np
T = {k: np.full(16, np.nan) for k in ("e_v", "e_3", "e_21", "e_22", "e_31", "e_4")}
for line in open(sys.argv[1]):
    if line.startswith("[gpkrec]"):
        t = line.split(); li = int(t[2])
        for k, v in zip(t[3::2], t[4::2]):
            T[k][li] = float(v)
np.savez(sys.argv[2], **T)
print("saved", sys.argv[2], {k: int(np.sum(~np.isnan(v))) for k, v in T.items()})
