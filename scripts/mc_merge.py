"""Merge Monte Carlo chunk files: python scripts/mc_merge.py OUT.npz IN1.npz IN2.npz ...  (sums add; N adds)."""
import sys, numpy as np
out = sys.argv[1]; acc = {}; N = 0
for f in sys.argv[2:]:
    z = np.load(f)
    for k in z.files:
        if k == "N": N += int(z["N"])
        elif k == "layers": acc["layers"] = z[k]
        else: acc[k] = acc.get(k, 0) + z[k]
np.savez(out, N=N, **acc); print(out, "N =", N, "keys", len(acc))
