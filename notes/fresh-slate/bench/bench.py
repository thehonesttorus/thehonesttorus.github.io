"""Shared benchmark sets for the fresh-slate design streams (numpy only; no flopscope, no parquet needed).

    import bench
    bench.list_sets()                       # names, width, depth, n_mlps, N
    S = bench.load_set("w128_d16")          # dict: seeds, width, depth, N, avg_variance (m,), means (m, L, n) float64,
                                            #       noise (m,) = avg_variance / N  (truth noise floor of the final-layer MSE)
    W = bench.weights(S, i)                 # (L, n, n) float32, x @ W convention, bit-identical to `whest dataset bake`
    for i, W in bench.iter_weights(S): ...

Weights are regenerated from the per-MLP seeds exactly as whestbench 0.16.1 does (seed protocol 3.0):
    weight_ss = SeedSequence(seed).spawn(3)[0];  rng = default_rng(weight_ss)
    W_l = (rng.standard_normal((n, n)) * sqrt(2/n)).astype(float32),  l = 1..L in order.
`python bench.py --verify DIR` checks bit-identity against a baked parquet dataset.
Truth: `all_layer_means` of the bake (stored float32 by whest, kept here as float64), `avg_variance` = mean per-neuron
variance of the final layer; the MC noise contribution to a final-layer MSE is avg_variance / N.
"""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SETS = os.path.join(HERE, "sets")


def list_sets():
    out = []
    for f in sorted(os.listdir(SETS)):
        if f.endswith(".json"):
            m = json.load(open(os.path.join(SETS, f)))
            out.append(dict(name=f[:-5], width=m["width"], depth=m["depth"], n_mlps=len(m["seeds"]),
                            N=m["n_samples"], max_noise=max(m["avg_variance"]) / m["n_samples"]))
    return out


def load_set(name):
    m = json.load(open(os.path.join(SETS, name + ".json")))
    z = np.load(os.path.join(SETS, name + ".npz"))
    m["means"] = z["means"].astype(np.float64)
    m["avg_variance"] = np.asarray(m["avg_variance"], dtype=np.float64)
    m["noise"] = m["avg_variance"] / float(m["n_samples"])
    m["name"] = name
    return m


def weights_from_seed(seed, width, depth):
    weight_ss = np.random.SeedSequence(int(seed)).spawn(3)[0]
    rng = np.random.default_rng(weight_ss)
    scale = float(np.sqrt(2.0 / width))
    return np.stack([(rng.standard_normal((width, width)) * scale).astype(np.float32) for _ in range(depth)])


def weights(S, i):
    return weights_from_seed(S["seeds"][i], S["width"], S["depth"])


def iter_weights(S):
    for i in range(len(S["seeds"])):
        yield i, weights(S, i)


def import_bake(bake_dir, name):
    """Convert a `whest dataset bake` directory into sets/<name>.json + .npz (small files only)."""
    import glob
    import pyarrow.parquet as pq
    meta = json.load(open(os.path.join(bake_dir, "metadata.json")))
    files = sorted(glob.glob(os.path.join(bake_dir, "data", "*.parquet")))
    rows = []
    for f in files:
        t = pq.read_table(f, columns=["mlp_id", "mlp_seed", "all_layer_means", "avg_variance"])
        rows += t.to_pylist()
    rows.sort(key=lambda r: r["mlp_id"])
    seeds_file = None
    # whest stores the estimator seed in mlp_seed; the input seeds come from the bake's seed list
    info = dict(width=meta["width"], depth=meta["depth"], n_samples=meta["n_samples"],
                seed_protocol=meta.get("seed_protocol"),
                avg_variance=[float(r["avg_variance"]) for r in rows],
                estimator_seeds=[int(r["mlp_seed"]) for r in rows])
    return info, np.stack([np.asarray(r["all_layer_means"], dtype=np.float64) for r in rows])


def save_set(name, info, means, input_seeds):
    os.makedirs(SETS, exist_ok=True)
    info = dict(info, seeds=[int(s) for s in input_seeds])
    json.dump(info, open(os.path.join(SETS, name + ".json"), "w"), indent=1)
    np.savez_compressed(os.path.join(SETS, name + ".npz"), means=means)


def verify(bake_dir, input_seeds):
    import glob
    import pyarrow.parquet as pq
    meta = json.load(open(os.path.join(bake_dir, "metadata.json")))
    files = sorted(glob.glob(os.path.join(bake_dir, "data", "*.parquet")))
    t = pq.read_table(files[0], columns=["mlp_id", "weights"])
    ok = True
    for k in range(min(t.num_rows, 2)):
        r = t.slice(k, 1).to_pylist()[0]
        Wb = np.asarray(r["weights"], dtype=np.float32)
        Wr = weights_from_seed(input_seeds[r["mlp_id"]], meta["width"], meta["depth"])
        same = Wb.shape == Wr.shape and np.array_equal(Wb, Wr)
        ok &= same
        print(f"mlp {r['mlp_id']}: bit-identical = {same}")
    return ok


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--verify", help="baked dataset dir")
    ap.add_argument("--import-bake", help="baked dataset dir")
    ap.add_argument("--seeds", help="JSON list of the input seeds given to `whest dataset bake --mlp-seeds`")
    ap.add_argument("--name")
    a = ap.parse_args()
    seeds = json.load(open(a.seeds)) if a.seeds else None
    if a.verify:
        print("OK" if verify(a.verify, seeds) else "MISMATCH")
    if a.import_bake:
        assert verify(a.import_bake, seeds), "weights do not regenerate bit-identically"
        info, means = import_bake(a.import_bake, a.name)
        save_set(a.name, info, means, seeds)
        print(a.name, means.shape, "max noise", max(info["avg_variance"]) / info["n_samples"])
    if not (a.verify or a.import_bake):
        for s in list_sets():
            print(s)
