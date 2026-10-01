"""Stage Q runner: python stageq.py <bake_dir> <variant> [variant ...]  -> JSON lines on stdout.
Variants: closure | hd:A=<int|inf>:diag=<all|star|startri>:second=<0|1>"""
import sys, json, time, glob
import numpy as np, pyarrow.parquet as pq
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from hd import closure, hd

def load(d):
    t = pq.read_table(glob.glob(d + '/data/*.parquet')[0]).to_pydict()
    meta = json.load(open(d + '/metadata.json'))
    for i in range(len(t['weights'])):
        W = np.array(t['weights'][i], dtype=np.float64)
        mu = np.array(t['all_layer_means'][i], dtype=np.float64)
        yield i, W, mu, t['avg_variance'][i] / meta['n_samples']

def parse(v):
    if v == 'closure':
        return lambda W: closure(W)
    kw = dict(x.split('=') for x in v.split(':')[1:])
    A = None if kw.get('A', 'inf') == 'inf' else int(kw['A'])
    return lambda W: hd(W, A=A, diagrams=kw.get('diag', 'all'), second=kw.get('second', '0') == '1')

if __name__ == '__main__':
    d, variants = sys.argv[1], sys.argv[2:]
    for i, W, mu, noise in load(d):
        for v in variants:
            t0 = time.time(); est = parse(v)(W)
            e = ((est - mu) ** 2).mean(1)
            print(json.dumps(dict(data=d, mlp=i, n=W.shape[1], variant=v, final=e[-1], noise=noise,
                                  final_minus_noise=e[-1] - noise, layers=list(e), sec=time.time() - t0)), flush=True)
