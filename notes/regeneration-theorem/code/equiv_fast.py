# equivalence of regen4_fast (BLAS-compiled) and regen4 (einsum) structure by structure, random small inputs
import re, math, numpy as np
from scipy.special import ndtr
def load(path):
    src = open(path).read()
    g = {"np": np, "math": math, "re": re, "ndtr": ndtr, "SQ2PI": math.sqrt(2 * math.pi), "EMAX": 4, "n": 9}
    exec(src[src.index("def he("):src.index("def r2(")], g)
    return g
A, B = load("regen4.py"), load("regen4_fast.py")
rng = np.random.default_rng(3); n = 9
mats = {k: rng.normal(size=(n, n)) for k in ("W", "Co", "D21", "K31", "K22")}
for k in ("Co", "D21", "K31", "K22"):
    np.fill_diagonal(mats[k], 0.0)
mats["Co"] = 0.5 * (mats["Co"] + mats["Co"].T); mats["K22"] = 0.5 * (mats["K22"] + mats["K22"].T)
F = {(p, d): rng.normal(size=n) for p in range(1, 5) for d in range(6)}
worst = 0.0; cnt = 0
for rows, od in (("iiij", False), ("iijj", False), ("iiii", True)):
    for st in A["structures"](rows, True):
        ya = A["evaluate"](rows, st, F, mats, od); yb = B["evaluate"](rows, st, F, mats, od)
        worst = max(worst, float(np.max(np.abs(ya - yb)) / (np.max(np.abs(ya)) + 1e-300))); cnt += 1
print(f"{cnt} structures, worst relative difference {worst:.2e}, einsum fallbacks {B['FALLBACK'][0]}")
