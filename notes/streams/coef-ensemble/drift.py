"""Drift of the ensemble coefficients away from the leg-partition values vs width: |c(n) - c_leg| per coefficient,
rms over layer bands, fitted as A n^-q (reads results/summary.json)."""
import json, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
S = json.load(open(os.path.join(HERE, "results", "summary.json")))
ws = sorted(int(w) for w in S)
LEG = np.array([1, 3, 3, 1, 1, 1.5, 1.5]); NAMES = ["B0", "B1", "B2", "B3", "B4", "B5", "B6"]
print("rms over layers of (c_ens(n) - c_leg), tensor-space table; fit |drift| ~ n^-q; value at n = 1024")
print("coef band  | " + " ".join(f"n={w:<5d}" for w in ws) + " |    q   drift(1024)")
for k in [2, 3, 5, 6, 4, 1]:
    for lo, hi in [(1, 3), (4, 9), (10, 14)]:
        d = [np.sqrt(np.mean((np.array(S[str(w)]["coef"]["tensor"])[lo:hi + 1, k] - LEG[k]) ** 2)) for w in ws]
        q, c = np.polyfit(np.log(ws), np.log(d), 1)
        print(f"{NAMES[k]} {lo:>2}-{hi:<2} | " + " ".join(f"{v:7.3f}" for v in d) + f" | {-q:5.2f}  {np.exp(c) * 1024 ** q:7.3f}")
