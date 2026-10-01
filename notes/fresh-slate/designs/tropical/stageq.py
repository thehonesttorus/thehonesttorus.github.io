"""Stage Q evaluator: run an estimator on baked datasets, report raw final-layer MSE minus truth noise.
Usage: python stageq.py <estimator: tct0|...> <dataset_dir> [<dataset_dir> ...]"""
import sys, json, time, importlib, os
import numpy as np, pyarrow.parquet as pq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
est_name = sys.argv[1]
mod, fn = (est_name.split(':') + [None])[:2]
est = getattr(importlib.import_module(mod), fn or mod)
for d in sys.argv[2:]:
    meta = json.load(open(os.path.join(d, 'metadata.json')))
    N = meta.get('n_samples') or meta.get('generation', {}).get('n_samples')
    t = pq.read_table([os.path.join(d, 'data', f) for f in os.listdir(os.path.join(d, 'data'))][0]).to_pydict()
    rows = []
    for i in range(len(t['weights'])):
        Ws = [np.asarray(w, dtype=np.float32) for w in t['weights'][i]]
        truth = np.asarray(t['all_layer_means'][i], dtype=np.float64)
        t0 = time.time(); pred = est(Ws); dt = time.time() - t0
        noise = t['avg_variance'][i] / N
        mse_layers = ((pred - truth) ** 2).mean(1)
        rows.append((mse_layers[-1] - noise, noise, mse_layers, dt))
        print(f"{d} mlp{i}: raw_final={mse_layers[-1]-noise:.3e} noise={noise:.1e} L1={mse_layers[0]:.1e} L4={mse_layers[3]:.1e} L8={mse_layers[7]:.1e} L12={mse_layers[11]:.1e} t={dt:.1f}s", flush=True)
    r = np.array([x[0] for x in rows])
    print(f"== {d} {est_name}: mean raw_final {r.mean():.3e} +- {r.std(ddof=1)/np.sqrt(len(r)):.1e} (n_mlps={len(r)}, N={N})", flush=True)
