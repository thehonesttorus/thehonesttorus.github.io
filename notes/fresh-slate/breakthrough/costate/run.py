"""Run co-state variants on a bench set; one JSON line per (variant, mlp) appended to results/<set>.jsonl."""
import sys, os, json, time, argparse, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "../../bench"))
import bench, costate

VARIANTS = {
    "gauss": dict(A=-1),
    "full": dict(),
    "A0": dict(A=0), "A1": dict(A=1), "A2": dict(A=2), "A3": dict(A=3), "A5": dict(A=5), "A7": dict(A=7),
    "A0slice": dict(A=0, old="slice"), "A1slice": dict(A=1, old="slice"), "A3slice": dict(A=3, old="slice"),
    "nocoinc": dict(coinc=False),
    "rk_half": dict(rank=lambda a: None if a < 1 else max(32, 1024 // (2 * a))),
}

ap = argparse.ArgumentParser()
ap.add_argument("--set", required=True); ap.add_argument("--variants", required=True)
ap.add_argument("--mlps", default=None)
a = ap.parse_args()
S = bench.load_set(a.set)
mlps = range(len(S["seeds"])) if a.mlps is None else [int(x) for x in a.mlps.split(",")]
os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
fn = os.path.join(HERE, "results", f"{a.set}.jsonl")
for i in mlps:
    W = bench.weights(S, i).astype(np.float64)
    truth = S["means"][i]
    for v in a.variants.split(","):
        kw = dict(VARIANTS[v])
        if v == "rk_half":
            n = W.shape[1]; kw["rank"] = (lambda nn: (lambda age: None if age < 1 else max(16, nn // (2 * age))))(n)
        rec = []
        t0 = time.time()
        pred = costate.predict(W, record=rec, **kw)
        dt = time.time() - t0
        lay = ((pred - truth) ** 2).mean(1)
        row = dict(set=a.set, variant=v, mlp=i, mse=float(lay[-1]), raw=float(lay[-1] - S["noise"][i]),
                   noise=float(S["noise"][i]), all_layer=float(lay.mean()), per_layer=[float(x) for x in lay],
                   nprod=rec[-1]["nprod"], sec=dt)
        with open(fn, "a") as f:
            f.write(json.dumps(row) + "\n")
        print(f"{a.set} {v:9s} mlp {i}: raw {row['raw']:.3e}  nprod {row['nprod']}  {dt:.0f}s", flush=True)
