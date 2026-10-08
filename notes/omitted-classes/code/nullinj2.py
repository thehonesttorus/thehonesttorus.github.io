# Part-by-part split of the V43 audit: which part of the null tuple fails to cancel, and where.
#   python nullinj2.py NET DIR
import sys, glob, os, re, numpy as np
net, d = int(sys.argv[1]), sys.argv[2]
R = {}
for f in glob.glob(f"{d}/ni_off{net}_L*_P*_+.npy"):
    m = re.search(r"_L(\d+)_(\w+?)_P(\d+)_\+\.npy$", f)
    fm = f.replace("_+.npy", "_-.npy")
    if os.path.exists(fm):
        R[(int(m.group(1)), m.group(2), int(m.group(3)))] = (np.load(f).astype(np.float64) - np.load(fm).astype(np.float64)) / 2
print(f"net {net}: output responses (last layer), as fractions of the covariance-part response r1; increments by part")
print("  layer gtype |  |r1|     |r3|/|r1| |r7|/|r1| |r15|/|r1| | kappa3 part |r3-r1|  kappa4 part |r7-r3|  feed |r15-r7|  (all / |r1|)"
      "  cos(r3-r1, -r1) cos(r7-r3, -r1)")
for L_, g in sorted({(k[0], k[1]) for k in R}):
    r = {P: R[(L_, g, P)][-1] for P in (1, 3, 7, 15) if (L_, g, P) in R}
    c = np.linalg.norm(r[1]); nr = lambda v: np.linalg.norm(v) / c
    cs = lambda a, b: float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-300))
    print(f"   {L_:2d}  {g:4s} | {c:.3e}  {nr(r[3]):.3f}     {nr(r[7]):.3f}     {nr(r[15]):.3f}    |   {nr(r[3] - r[1]):.3f}"
          f"             {nr(r[7] - r[3]):.3f}             {nr(r[15] - r[7]):.3f}           {cs(r[3] - r[1], -r[1]):+.3f}          {cs(r[7] - r[3], -r[1]):+.3f}")
