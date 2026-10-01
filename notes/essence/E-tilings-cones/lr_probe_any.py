"""Copy of region's lr_probe.py (mean and off probes only) for any bench set, to test the wedge-calculus prediction
(Proposition E1) away from 1024 x 16.  usage: python lr_probe_any.py SET MLP"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, "..", "..", "fresh-slate", "breakthrough", "region")
sys.path.insert(0, os.path.join(HERE, "..", "..", "fresh-slate", "bench")); sys.path.insert(0, R)
import bench, gclose as g

def main(name, i):
    S = bench.load_set(name)
    W = bench.weights(S, i).astype(np.float64)
    L, n, _ = W.shape
    base, states = g.run(W, "exact")
    rng = np.random.default_rng(77 + i)
    res = {"mean": [None] * L, "off": [None] * L, "n": n, "L": L}
    for l in range(L - 1):
        m0, C0 = states[l]
        for p in ("mean", "off"):
            if p == "mean":
                dm = rng.standard_normal(n); dC = None; amp = 1e-3
            else:
                X = rng.standard_normal((n, n)); X = (X + X.T) / np.sqrt(2); np.fill_diagonal(X, 0)
                dm = None; dC = X / np.linalg.norm(X); amp = 1e-2
            outs = []
            for sgn in (+1, -1):
                m1 = m0 + (sgn * amp * dm if dm is not None else 0)
                C1 = C0 + (sgn * amp * dC if dC is not None else 0)
                o, _ = g.run(W, "exact", start=(l + 1, m1, C1)); outs.append(o[-1])
            d = (outs[0] - outs[1]) / (2 * amp)
            norm2 = float(np.mean(dm ** 2)) if p == "mean" else 1.0
            res[p][l] = float(np.mean(d ** 2) / norm2)
        print(f"{name} mlp {i} l {l}: Kmean {res['mean'][l]:.3e} Koff {res['off'][l]:.3e}", flush=True)
    json.dump(res, open(os.path.join(HERE, "results", f"lr_{name}_mlp{i}.json"), "w"), indent=1)

if __name__ == "__main__":
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    main(sys.argv[1], int(sys.argv[2]))
