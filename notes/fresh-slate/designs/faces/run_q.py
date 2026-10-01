"""Stage Q driver: bench sets via eval_q, plus local bakes (same format) for sets not yet on the bench."""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "bench")); sys.path.insert(0, HERE)
import bench, eval_q, fbt

def local_set(path, name):
    import pyarrow.parquet as pq, glob
    rows = []
    for f in sorted(glob.glob(os.path.join(path, "data", "*.parquet"))):
        rows += pq.read_table(f, columns=["mlp_seed", "all_layer_means", "avg_variance"]).to_pylist()
    meta = json.load(open(os.path.join(path, "metadata.json")))
    S = dict(seeds=[r["mlp_seed"] for r in rows], width=meta["width"], depth=meta["depth"], n_samples=meta["n_samples"],
             means=np.array([r["all_layer_means"] for r in rows], dtype=np.float64),
             avg_variance=np.array([r["avg_variance"] for r in rows]), name=name)
    S["noise"] = S["avg_variance"] / S["n_samples"]
    return S

_orig = bench.load_set
def load_set(name):
    p = os.path.join("/root/bench", name.split("_")[0])
    if name.startswith("my") :
        return local_set(os.path.join("/root/bench", "w" + name[2:].split("_")[0]), name)
    return _orig(name)
bench.load_set = load_set

if __name__ == "__main__":
    funcs = sys.argv[1].split(","); sets = sys.argv[2].split(","); units = float(sys.argv[3]) if len(sys.argv) > 3 else None
    allres = {}
    for fn in funcs:
        print(f"## {fn}", flush=True)
        res = eval_q.evaluate(getattr(fbt, fn), sets, units_1024=units)
        allres[fn] = res
    out = os.path.join(HERE, "results", f"q_{'_'.join(funcs)}_{'_'.join(sets)}.json")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    json.dump(allres, open(out, "w"), indent=1)
