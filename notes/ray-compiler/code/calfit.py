# How much of the chain's error can any per-layer amplitude correction of its carried statistics remove?
# (ray-compiler note section 8.) Linear responses R_j = (out_cj - out_b)[-1] / d_j of the final means to rescaling
# statistic X_j at layer l_j by 1 + d_j (V47_CAL), then the output-metric projection b = argmin sum_train |e + R b|^2
# + ridge |b|^2 (ridge chosen inside the training set), and the predicted held-out MSE change.
#   python calfit.py RESPDIR BASEDIR DIRS_JSON [--train 0-49] [--val 50-99]
# Also reports the projection per statistic family (only that family's directions) and per layer band.
import sys, os, json, numpy as np


def nets(s):
    a, _, b = s.partition("-"); return list(range(int(a), int(b) + 1)) if b else [int(a)]


rd, bd, dj = sys.argv[1], sys.argv[2], sys.argv[3]
tr = nets(sys.argv[sys.argv.index("--train") + 1]) if "--train" in sys.argv else list(range(50))
va = nets(sys.argv[sys.argv.index("--val") + 1]) if "--val" in sys.argv else list(range(50, 100))
dirs = json.load(open(dj)); K = len(dirs)
off = os.environ.get("OFFICIAL", "../official")


def load(n):
    m = np.load(f"{off}/truth_off{n}.npz")["m"].astype(np.float64)
    b = np.load(f"{bd}/out_b_{n}.npy")[-1]
    R = np.stack([(np.load(f"{rd}/out_c{j}_{n}.npy")[-1] - b) / dirs[j][2] for j in range(K)], axis=1)
    return b - m[-1], R


data = {n: load(n) for n in tr + va if all(os.path.exists(f"{rd}/out_c{j}_{n}.npy") for j in range(K))}
tr = [n for n in tr if n in data]; va = [n for n in va if n in data]
print(f"{K} directions; train {len(tr)} nets, held-out {len(va)} nets")


def fit(idx, ridges=(0.0, 1e-4, 1e-3, 1e-2, 1e-1, 1.0)):
    idx = np.asarray(idx)
    def gram(ns):
        G = sum(data[n][1][:, idx].T @ data[n][1][:, idx] for n in ns); g = sum(data[n][1][:, idx].T @ data[n][0] for n in ns)
        return G, g
    h1, h2 = tr[:len(tr) // 2], tr[len(tr) // 2:]
    G1, g1 = gram(h1)
    def sc(rg):
        bb = -np.linalg.solve(G1 + rg * np.trace(G1) / len(idx) * np.eye(len(idx)), g1)
        return np.mean([np.mean((data[n][0] + data[n][1][:, idx] @ bb) ** 2) for n in h2])
    rg = min(ridges, key=sc)
    G, g = gram(tr)
    b = -np.linalg.solve(G + rg * np.trace(G) / len(idx) * np.eye(len(idx)), g)
    def ch(ns):
        m0 = np.array([np.mean(data[n][0] ** 2) for n in ns]); m1 = np.array([np.mean((data[n][0] + data[n][1][:, idx] @ b) ** 2) for n in ns])
        return 100 * (m1.mean() / m0.mean() - 1), 100 * (m1 / m0 - 1).std() / np.sqrt(len(ns)), int((m1 < m0).sum())
    return b, rg, ch(tr), ch(va)


b, rg, ctr, cva = fit(range(K))
print(f"ALL {K} directions: ridge {rg:g}; predicted train {ctr[0]:+.2f}%, held-out {cva[0]:+.2f}% +- {cva[1]:.2f} (better on {cva[2]}/{len(va)})")
print("  coefficients (statistic:layer:factor): " + " ".join(f"{x}:{l}:{1 + bb * 1.0:.4f}" if x in ('var', 'coff') else f"{x}:{l}:{1 + bb:.3f}"
                                                     for (x, l, _), bb in zip(dirs, b)))
for fam in ("var", "coff", "D3", "D21", "g4", "k22", "k31"):
    idx = [j for j, d in enumerate(dirs) if d[0] == fam]
    b2, rg2, c2, v2 = fit(idx)
    print(f"  {fam:5s} alone ({len(idx):2d} dirs): ridge {rg2:g}; train {c2[0]:+.2f}%, held-out {v2[0]:+.2f}% +- {v2[1]:.2f}")
for lo, hi in ((0, 4), (5, 9), (10, 15)):
    idx = [j for j, d in enumerate(dirs) if lo <= d[1] <= hi]
    b2, rg2, c2, v2 = fit(idx)
    print(f"  layers {lo:2d}-{hi:2d} ({len(idx):2d} dirs): train {c2[0]:+.2f}%, held-out {v2[0]:+.2f}% +- {v2[1]:.2f}")
json.dump({"dirs": dirs, "coef": [float(x) for x in b], "ridge": rg}, open(os.environ.get("CALFIT_OUT", "calfit_coef.json"), "w"))
