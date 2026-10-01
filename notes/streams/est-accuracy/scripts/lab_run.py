"""In-process runner for estimator variants on a baked dataset (one estimator instance across MLPs,
as in the suite).  Saves every layer's prediction so variants can be compared paired, per layer.

  python lab_run.py --estimator EST.py --dataset DIR --tag NAME [--env K=V ...] [--mlps 0,1,2]

Writes results/<tag>.json (per-MLP final MSE, MSE - noise, FLOPs/B, residual, wall) and
results/preds/<tag>.npz (preds (m, L, n) float32).  Weights are cached as .npy next to the dataset.
"""
import argparse, glob, json, os, sys, time

ap = argparse.ArgumentParser()
ap.add_argument("--estimator", required=True)
ap.add_argument("--dataset", required=True)
ap.add_argument("--tag", required=True)
ap.add_argument("--env", nargs="*", default=[])
ap.add_argument("--mlps", default="")
ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results"))
a = ap.parse_args()
for kv in a.env:
    k, v = kv.split("=", 1)
    os.environ[k] = v

import numpy as np
import importlib.util
import flopscope as flops
import flopscope.numpy as fnp
from whestbench import MLP, SetupContext


def load_dataset(d):
    cache = os.path.join(d, "npy")
    meta = json.load(open(os.path.join(d, "metadata.json")))
    if not os.path.exists(os.path.join(cache, "truth.npy")):
        import pyarrow.parquet as pq
        os.makedirs(cache, exist_ok=True)
        files = sorted(glob.glob(os.path.join(d, "data", "*.parquet")))
        tab = pq.read_table(files[0])
        truths, avgv, seeds = [], [], []
        for i in range(tab.num_rows):
            r = tab.slice(i, 1).to_pylist()[0]
            np.save(os.path.join(cache, f"w{i}.npy"), np.asarray(r["weights"], dtype=np.float32))
            truths.append(np.asarray(r["all_layer_means"], dtype=np.float32))
            avgv.append(float(r["avg_variance"]))
            seeds.append(int(r["mlp_seed"]))
        np.save(os.path.join(cache, "truth.npy"), np.stack(truths))
        json.dump(dict(avg_variance=avgv, seeds=seeds), open(os.path.join(cache, "info.json"), "w"))
    truth = np.load(os.path.join(cache, "truth.npy"))
    info = json.load(open(os.path.join(cache, "info.json")))
    return meta, truth, info, cache


meta, truth, info, cache = load_dataset(a.dataset)
N = float(meta["n_samples"])
idx = [int(x) for x in a.mlps.split(",")] if a.mlps else list(range(truth.shape[0]))

spec = importlib.util.spec_from_file_location("estimator", a.estimator)
mod = importlib.util.module_from_spec(spec)
sys.path.insert(0, os.path.dirname(os.path.abspath(a.estimator)))
spec.loader.exec_module(mod)
est = mod.Estimator()
with flops.BudgetContext(flop_budget=10 ** 15, quiet=True):
    est.setup(SetupContext(width=int(meta["width"]), depth=int(meta["depth"]), flop_budget=2 ** 41,
                           api_version="1.0", submission_dir=os.path.dirname(os.path.abspath(a.estimator)), seed=0))
B = 2 ** 41
rows, preds = [], []
for i in idx:
    w = np.load(os.path.join(cache, f"w{i}.npy"))
    with flops.BudgetContext(flop_budget=10 ** 15, quiet=True):
        weights = [fnp.asarray(w[l]) for l in range(w.shape[0])]
    mlp = MLP(width=w.shape[1], depth=w.shape[0], weights=weights, seed=int(info["seeds"][i]) % (2 ** 63))
    ctx = flops.BudgetContext(flop_budget=B, quiet=True)
    rec = dict(mlp=i)
    t0 = time.time()
    try:
        with ctx:
            pred = est.predict(mlp, B)
            fl = ctx.flops_used
        p = np.asarray(pred, dtype=np.float32)
        noise = info["avg_variance"][i] / N
        mse = float(np.mean((p[-1] - truth[i, -1]) ** 2))
        rec.update(ok=True, cb=fl / B, mse=mse, noise=noise, adj=mse - noise,
                   layer_mse=[float(np.mean((p[l] - truth[i, l]) ** 2)) for l in range(p.shape[0])],
                   finite=bool(np.isfinite(p).all()))
        preds.append(p)
    except Exception as e:  # noqa
        rec.update(ok=False, error=f"{type(e).__name__}: {e}"[:400])
        preds.append(np.full(truth[i].shape, np.nan, np.float32))
    rec.update(wall=time.time() - t0, residual=ctx.residual_wall_time_s)
    if getattr(mod, "EA_DUMP", ""):
        dd = {}
        for d in mod.EA_DUMPS:
            for k, v in d.items():
                if k != "li" and v is not None:
                    dd[f"{k}_{d['li']}"] = v
        os.makedirs(mod.EA_DUMP, exist_ok=True)
        np.savez(os.path.join(mod.EA_DUMP, f"{a.tag}_mlp{i}.npz"), **dd)
        mod.EA_DUMPS.clear()
    rows.append(rec)
    print(json.dumps({k: v for k, v in rec.items() if k != "layer_mse"}), flush=True)
ok = [r for r in rows if r.get("ok")]
summ = dict(tag=a.tag, env=a.env, estimator=a.estimator, dataset=a.dataset, n_samples=N, mlps=idx,
            mean_mse=float(np.mean([r["mse"] for r in ok])) if ok else None,
            mean_adj=float(np.mean([r["adj"] for r in ok])) if ok else None,
            cb=float(np.max([r["cb"] for r in ok])) if ok else None,
            residual_max=float(np.max([r["residual"] for r in rows])), per_mlp=rows)
os.makedirs(os.path.join(a.out, "preds"), exist_ok=True)
json.dump(summ, open(os.path.join(a.out, f"{a.tag}.json"), "w"), indent=1)
np.savez_compressed(os.path.join(a.out, "preds", f"{a.tag}.npz"), preds=np.stack(preds), mlps=np.array(idx))
print("SUMMARY", a.tag, summ["mean_mse"], summ["mean_adj"], summ["cb"], summ["residual_max"])
