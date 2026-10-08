# Output-metric refit of the lam table (ray-compiler note section 8).
#   python lamfit.py DIR [--train 0-49] [--val 50-99] [--ridge R]
# DIR holds out_b_N.npy (production, all-layer means) and out_pL_N.npy (LAM[L] x 1.1, L = 0..14) for every network.
# The chain's final-layer error is e = out_b[-1] - truth[-1]. A small change of the multipliers, m_L = 1 + b_L, moves the
# output by sum_L b_L R_L with R_L = (out_pL - out_b)[-1] / 0.1 (linear response along the chain's own trajectory).
# The output-metric projection is b = argmin sum_train |e + R b|^2 + ridge |b|^2, so the predicted held-out MSE is
# mean_val |e + R b|^2 / n: the projection of the one-step defects onto the lam direction, measured where the output
# reads it.
import sys, os, numpy as np


def nets(s):
    a, _, b = s.partition("-"); return list(range(int(a), int(b) + 1)) if b else [int(a)]


d = sys.argv[1]
tr = nets(sys.argv[sys.argv.index("--train") + 1]) if "--train" in sys.argv else list(range(50))
va = nets(sys.argv[sys.argv.index("--val") + 1]) if "--val" in sys.argv else list(range(50, 100))
ridges = [float(sys.argv[sys.argv.index("--ridge") + 1])] if "--ridge" in sys.argv else [0.0, 1e-3, 1e-2, 1e-1]
off = os.environ.get("OFFICIAL", "../official")
K = 15


def load(n):
    m = np.load(f"{off}/truth_off{n}.npz")["m"].astype(np.float64)
    b = np.load(f"{d}/out_b_{n}.npy")
    e = b[-1] - m[-1]
    R = np.stack([(np.load(f"{d}/out_p{l}_{n}.npy")[-1] - b[-1]) / 0.1 for l in range(K)], axis=1)   # (n, K)
    return e, R


data = {n: load(n) for n in tr + va if os.path.exists(f"{d}/out_b_{n}.npy")}
tr = [n for n in tr if n in data]; va = [n for n in va if n in data]
G = sum(data[n][1].T @ data[n][1] for n in tr); g = sum(data[n][1].T @ data[n][0] for n in tr)
scale = np.trace(G) / K
print(f"train {len(tr)} nets, validation {len(va)} nets; response norms per layer (rms over train): "
      + " ".join(f"{np.sqrt(G[l, l] / len(tr) / 1024):.1e}" for l in range(K)))
# ridge chosen inside the training set (first half fits, second half scores), never on the held-out networks
h1, h2 = tr[:len(tr) // 2], tr[len(tr) // 2:]
G1 = sum(data[n][1].T @ data[n][1] for n in h1); g1 = sum(data[n][1].T @ data[n][0] for n in h1)
def _score(rg):
    bb = -np.linalg.solve(G1 + rg * np.trace(G1) / K * np.eye(K), g1)
    return np.mean([np.mean((data[n][0] + data[n][1] @ bb) ** 2) for n in h2])
best = min(ridges, key=_score)
print("ridge chosen inside the training set: " + " ".join(f"{rg:g}:{_score(rg):.4e}" for rg in ridges) + f" -> {best:g}")
for rg in ridges:
    b = -np.linalg.solve(G + rg * scale * np.eye(K), g)
    def mse(ns, bb):
        return np.array([np.mean((data[n][0] + data[n][1] @ bb) ** 2) for n in ns])
    m0t, m1t = mse(tr, 0 * b), mse(tr, b); m0v, m1v = mse(va, 0 * b), mse(va, b)
    rv = m1v / m0v - 1
    print(f"ridge {rg:g}{' (chosen)' if rg == best else ''}: multipliers " + " ".join(f"{1 + x:.3f}" for x in b))
    print(f"   predicted (first order): train {100 * (m1t.mean() / m0t.mean() - 1):+.2f}%, "
          f"held-out {100 * (m1v.mean() / m0v.mean() - 1):+.2f}% (per net {100 * rv.mean():+.2f} +- {100 * rv.std() / np.sqrt(len(rv)):.2f}, "
          f"better on {int((rv < 0).sum())}/{len(rv)})")
