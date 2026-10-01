"""Stage Q evaluation against a whest bake: raw MSE per layer minus truth noise."""
import sys, time, numpy as np, pyarrow.parquet as pq, glob
sys.path.insert(0, '.')
import bethe

def load(path):
    rows = []
    for f in sorted(glob.glob(path + '/data/*.parquet')):
        d = pq.read_table(f).to_pydict()
        for i in range(len(d['mlp_id'])):
            rows.append(dict(W=np.array(d['weights'][i]), mean=np.array(d['all_layer_means'][i]),
                             avgvar=d['avg_variance'][i]))
    return rows

def main(path, N, methods):
    rows = load(path)
    res = {k: [] for k in methods}
    for r in rows:
        noise = r['avgvar'] / N   # final-layer truth noise (approx, all layers similar order)
        for k in methods:
            t = time.time()
            e = getattr(bethe, k[0])(r['W'], **k[1]) if isinstance(k, tuple) else getattr(bethe, k)(r['W'])
            mse = ((e - r['mean']) ** 2).mean(1)
            res[k].append(mse[-1] - noise)
    for k, v in res.items():
        print(f"{str(k):40s} final MSE-noise: mean {np.mean(v):.3e}  per-MLP {' '.join(f'{x:.2e}' for x in v)}")

if __name__ == '__main__':
    path, N = sys.argv[1], float(sys.argv[2])
    ms = ['estimate_tree', 'estimate_gauss', ('estimate_fact', dict(old=0)), ('estimate_fact', dict(old=1)), ('estimate_full', {})]
    ms = [m if isinstance(m, str) else (m[0], m[1]) for m in ms]
    rows = load(path)
    for r in rows:
        noise = r['avgvar'] / N
        line = []
        for m in ms:
            f = getattr(bethe, m) if isinstance(m, str) else (lambda W, m=m: getattr(bethe, m[0])(W, **m[1]))
            e = f(r['W'])
            mse = ((e - r['mean']) ** 2).mean(1) - noise
            line.append(f"{mse[-1]:.2e}")
        print('noise %.1e | tree gauss fact0 fact1 full:' % noise, ' '.join(line), flush=True)
