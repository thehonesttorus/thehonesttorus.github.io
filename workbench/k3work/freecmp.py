# Free-running chain dumps against Monte Carlo truth, per layer and channel (note XXXVI).
#   python freecmp.py NET MCPREFIX DUMP1 DUMP2
import sys, numpy as np
net, pre, d1, d2 = int(sys.argv[1]), sys.argv[2], sys.argv[3], sys.argv[4]
F = np.load(f"{pre}_full.npz"); n = F["mu"].shape[1]; off = ~np.eye(n, dtype=bool)
A, B = np.load(d1), np.load(d2)
def get(ch, k, l):
    key = f"{k}_{l}"
    return None if key not in ch.files else ch[key].astype(np.float64)
tr = {"var": ("var", None), "D3": ("k3", None), "D21": ("D21", off), "g4row": ("k4", None), "wk4m": ("K22", off), "C_off": ("cov", off)}
print(f"net {net}: free-running relative error vs truth, {d1} | {d2}  (and corr of the two error vectors)")
for l in range(1, 16):
    row = []
    for k, (tk, msk) in tr.items():
        a, b = get(A, k, l), get(B, k, l)
        if a is None or b is None:
            continue
        t = F[tk][l].astype(np.float64)
        if msk is not None:
            a, b, t = a[msk], b[msk], t[msk]
        ea, eb = a - t, b - t
        row.append(f"{k} {np.linalg.norm(ea) / np.linalg.norm(t):.4f}|{np.linalg.norm(eb) / np.linalg.norm(t):.4f} ({np.corrcoef(ea, eb)[0, 1]:+.2f})")
    print(f"layer {l:2d}: " + "  ".join(row), flush=True)
