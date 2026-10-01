"""Aggregate results/raw/*.json into a table: final-layer MSE per variant (mean over MLPs, raw and minus the truth
noise floor avg_variance / N), plus per-layer MSE at a few depths.
    python summarize.py [--dataset /root/work/truth128] [--mlps 0-7] > results/summary.md
"""
import argparse, glob, json, os
import numpy as np
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", default="results/raw"); ap.add_argument("--dataset", default="/root/work/truth128")
    ap.add_argument("--mlps", default=None, help="comma list; default = all MLPs present for the variant")
    a = ap.parse_args()
    floor = {}
    try:
        import pyarrow.parquet as pq
        meta = json.load(open(os.path.join(a.dataset, "metadata.json")))
        N = float(meta["n_samples"])
        t = pq.read_table(glob.glob(os.path.join(a.dataset, "data", "*.parquet"))[0], columns=["avg_variance"])
        for i, v in enumerate(t.column("avg_variance").to_pylist()):
            floor[i] = v / N
    except Exception:
        pass
    want = None if a.mlps is None else set(int(x) for x in a.mlps.split(","))
    res = defaultdict(dict)
    for f in glob.glob(os.path.join(a.raw, "*.json")):
        d = json.load(open(f))
        if want is None or d["mlp"] in want:
            res[d["variant"]][d["mlp"]] = d
    print(f"| variant | MLPs | final MSE (mean) | minus floor | median | layer 1 | layer 4 | layer 8 | layer 12 |")
    print("|---|---|---|---|---|---|---|---|---|")
    rows = []
    for v, dd in res.items():
        mlps = sorted(dd)
        fin = np.array([dd[i]["final_mse"] for i in mlps])
        fl = np.array([floor.get(i, 0.0) for i in mlps])
        pl = np.array([dd[i]["per_layer_mse"] for i in mlps]).mean(0)
        rows.append((fin.mean(), v, mlps, fin, fl, pl))
    for m, v, mlps, fin, fl, pl in sorted(rows, key=lambda r: r[0]):
        ids = ",".join(map(str, mlps))
        print(f"| {v} | {ids} | {fin.mean():.3e} | {(fin - fl).mean():.3e} | {np.median(fin):.3e} | {pl[1]:.2e} | {pl[4]:.2e} | {pl[8]:.2e} | {pl[12]:.2e} |")
    if floor:
        print(f"\ntruth noise floor avg_variance/N: mean {np.mean(list(floor.values())):.2e}")


if __name__ == "__main__":
    main()
