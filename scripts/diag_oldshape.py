"""Are the old kappa3 sources' slice contributions shape-similar to the young ones? Dense two-tier chain (no merging);
per layer: corr and least-squares scale of the old (age >= w) (2,1)-slice contribution against the age-1 contribution
and against the sum of ages < w, entrywise over the n^2 off-diagonal entries and over the diagonal (3,) slices.
  python scripts/diag_oldshape.py DATA NET [w=2]"""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from whest import kprop3c
from whest.kprop3c import kprop3c_chain, slices_from_legs
D, net = sys.argv[1], int(sys.argv[2]); w = int(sys.argv[3].split("=")[1]) if len(sys.argv) > 3 else 2
W = np.load(f"{D}/W_off{net}.npy").astype(np.float64)
per_age = {}
orig = kprop3c.nonlin_step
def hooked(st, o, rng, record=None):
    l = len(per_age); per_age[l] = {age: slices_from_legs(legs) for legs, age in st.young}
    return orig(st, o, rng, record)
kprop3c.nonlin_step = hooked
kprop3c_chain(W, dict(window=99, k=8))
def stats(X, Y):
    x, y = X.ravel(), Y.ravel(); a = np.dot(x, y) / np.dot(x, x); return np.corrcoef(x, y)[0, 1], a, np.sqrt(np.mean((y - a * x) ** 2) / np.mean(y ** 2))
print("layer | old(age>=w) vs age1: corr, scale, resid | old vs young(<w): corr, scale, resid | same for the (3,) diagonal | |old|/|young| (21)")
for l in sorted(per_age):
    d = per_age[l]
    if not any(a >= w for a in d): continue
    old21 = sum(d[a][0] for a in d if a >= w); old3 = sum(d[a][1] for a in d if a >= w)
    y21 = sum(d[a][0] for a in d if a < w); y3 = sum(d[a][1] for a in d if a < w)
    c1, a1, r1 = stats(d[1][0], old21); c2, a2, r2 = stats(y21, old21); c3, a3, r3 = stats(y3, old3)
    print(f"  {l:2d} | {c1:+.3f} {a1:+.3f} {r1:.3f} | {c2:+.3f} {a2:+.3f} {r2:.3f} | {c3:+.3f} {a3:+.3f} {r3:.3f} | {np.linalg.norm(old21)/np.linalg.norm(y21):.3f}", flush=True)
