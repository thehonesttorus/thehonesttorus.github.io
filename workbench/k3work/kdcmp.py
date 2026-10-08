# One-step kappa4 slices of the chain (own values in an all-oracle dump) against Monte Carlo truth.
#   python kdcmp.py NET MCPREFIX DUMP1 [DUMP2 ...]
import sys, numpy as np
net, pre, dumps = int(sys.argv[1]), sys.argv[2], sys.argv[3:]
F = np.load(f"{pre}_full.npz"); n = F["mu"].shape[1]
rng = np.random.default_rng(1); ia = rng.integers(0, n, 12000); ib = rng.integers(0, n, 12000); k = ia != ib; ia, ib = ia[k], ib[k]
off = ~np.eye(n, dtype=bool)
chs = [np.load(d) for d in dumps]
print(f"net {net}: one-step relative error diag | (2,2) | (3,1) of " + " vs ".join(dumps))
for l in range(1, 15):
    t = (F["k4"][l + 1].astype(np.float64), F["K22"][l + 1].astype(np.float64)[ia, ib], F["K31"][l + 1].astype(np.float64)[off])
    row = []
    for ch in chs:
        if any(f"{k}own_{l + 1}" not in ch.files for k in ("g4row", "wk4m", "wk431")):
            row.append("-"); continue
        o = (ch[f"g4rowown_{l + 1}"].astype(np.float64), ch[f"wk4mown_{l + 1}"].astype(np.float64)[ia, ib],
             ch[f"wk431own_{l + 1}"].astype(np.float64).T[off])
        row.append(" ".join(f"{np.linalg.norm(a - b) / np.linalg.norm(b):.3f}" for a, b in zip(o, t)))
    print(f"layer {l:2d}->{l + 1:2d}: " + " | ".join(row), flush=True)
