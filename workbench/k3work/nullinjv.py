# The V43 null audit under mechanism switches: which part of the chain carries the anomaly.
#   python nullinjv.py NET LABEL=DIR [LABEL=DIR ...]
import sys, os, numpy as np
net = int(sys.argv[1]); runs = [a.split("=", 1) for a in sys.argv[2:]]
def resp(d, L_, g, P):
    f = f"{d}/ni_off{net}_L{L_}_{g}_P{P}_"
    if not (os.path.exists(f + "+.npy") and os.path.exists(f + "-.npy")):
        return None
    return (np.load(f + "+.npy").astype(np.float64) - np.load(f + "-.npy").astype(np.float64)) / 2
print(f"net {net}: |r_full| / |r_cov| at the output (exact 0), and |r_cov| (the covariance part's own response)")
print("  " + " ".join(f"{'L%d %s' % (L_, g):>22s}" for L_ in (7, 11) for g in ("var", "rnd")))
for lab, d in runs:
    cells = []
    for L_ in (7, 11):
        for g in ("var", "rnd"):
            r1, r15 = resp(d, L_, g, 1), resp(d, L_, g, 15)
            cells.append("        n/a           " if r1 is None or r15 is None else
                         f"{np.linalg.norm(r15[-1]) / np.linalg.norm(r1[-1]):6.3f} ({np.linalg.norm(r1[-1]):.2e})  ")
    print(f"  {lab:28s} " + " ".join(cells))
