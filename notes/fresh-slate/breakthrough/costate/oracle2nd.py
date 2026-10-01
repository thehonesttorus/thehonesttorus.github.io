"""Second-order co-state oracle (w128 bench nets): true MC cumulant slices injected on the own (m, C) chain.
Question: how much does true joint kappa_4 gain over true kappa_3 (first order), and is its readout-relevant
part concentrated (rank-r off-diagonal + exact diagonal)?  Atlas code: heisenberg oracle3 (two passes)."""
import sys, os, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "../../bench")); sys.path.insert(0, os.path.join(HERE, "../../designs/heisenberg"))
import bench
from hd import Gauss, inject, inject_full2


def passes(Ws, N, seed, mean=None, chunk=50000):
    L, n, _ = Ws.shape; rng = np.random.default_rng(seed); done = 0
    acc = dict(m=np.zeros((L, n)), a=np.zeros((L, n)))
    if mean is not None:
        for k in ["c2", "s21", "s22", "s31"]: acc[k] = np.zeros((L, n, n))
        acc["d3"] = np.zeros((L, n)); acc["d4"] = np.zeros((L, n))
    while done < N:
        b = min(chunk, N - done); a = rng.standard_normal((b, n))
        for l in range(L):
            z = a @ Ws[l]; a = np.maximum(z, 0); acc["a"][l] += a.sum(0); acc["m"][l] += z.sum(0)
            if mean is not None:
                x = z - mean[l]; x2 = x * x
                acc["c2"][l] += x.T @ x; acc["s21"][l] += x2.T @ x; acc["s22"][l] += x2.T @ x2; acc["s31"][l] += (x2 * x).T @ x
                acc["d3"][l] += (x2 * x).sum(0); acc["d4"][l] += (x2 * x2).sum(0)
        done += b
    return {k: v / N for k, v in acc.items()}


def lowrank_off(X, r):
    if r is None: return X
    d = np.diag(X).copy(); O = X - np.diag(d)
    if r == 0: return np.diag(d)
    U, s, Vt = np.linalg.svd(O); return (U[:, :r] * s[:r]) @ Vt[:r] + np.diag(d)


def run(Ws, st, mode, r=None):
    L, n, _ = Ws.shape; m = np.zeros(n); C = Ws[0].T @ Ws[0]; res = []
    for l in range(L):
        G = Gauss(m, C, K=21)
        c2 = st["c2"][l]; v = np.diag(c2)
        D = st["d3"][l]; S = st["s21"][l]; K4d = st["d4"][l] - 3 * v * v
        K22 = st["s22"][l] - np.outer(v, v) - 2 * c2 * c2; K31 = st["s31"][l] - 3 * v[:, None] * c2
        if l == 0: dEa, dC = np.zeros(n), np.zeros((n, n))
        elif mode == "first": dEa, dC = inject(G, D, S)
        elif mode == "full2": dEa, dC = inject_full2(G, D, S, K4d, lowrank_off(K22, r), lowrank_off(K31, r))
        elif mode == "k4diag": dEa, dC = inject_full2(G, D, S, K4d, np.diag(np.diag(K22)), np.diag(np.diag(K31)))
        Ea = G.Ea + dEa; res.append(Ea)
        if l + 1 < L: Ca = G.cov_a() + dC; m = Ea @ Ws[l + 1]; C = Ws[l + 1].T @ Ca @ Ws[l + 1]
    return np.array(res)


if __name__ == "__main__":
    name, N = sys.argv[1], int(float(sys.argv[2])); mlps = [int(x) for x in sys.argv[3].split(",")]
    S = bench.load_set(name)
    out = os.path.join(HERE, "results_live", f"oracle2nd_{name}.jsonl"); os.makedirs(os.path.dirname(out), exist_ok=True)
    for i in mlps:
        Ws = bench.weights(S, i).astype(np.float64)
        p1 = passes(Ws, N, 1000 + i); st = passes(Ws, N, 2000 + i, mean=p1["m"])
        for mode, r in [("first", None), ("full2", None), ("full2", 1), ("full2", 8), ("full2", 32), ("k4diag", None)]:
            pred = run(Ws, st, mode, r)
            raw = float(((pred[-1] - S["means"][i][-1]) ** 2).mean() - S["noise"][i])
            row = dict(set=name, mlp=i, mode=mode, r=r, N=N, raw=raw)
            open(out, "a").write(json.dumps(row) + "\n")
            print(name, i, mode, r, f"raw {raw:.3e}", flush=True)
