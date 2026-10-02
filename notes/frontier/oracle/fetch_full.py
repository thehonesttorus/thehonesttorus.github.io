"""Fetch N=1e9 truth for the public `full` split (needs huggingface.co + *.hf.co allowed by the network policy).
1. weights + all_layer_means: one parquet shard of aicrowd/arc-whestbench-public-2026 (rev v2-phase2, ~15 MLPs, 1.1 GB)
2. per-layer pre-activation marginals of the same networks: keenanpepper/arc-whestbench-p2-full1000-N1e9, range-read
   members pre_mean, pre_m2, pre_m3, pre_m4, gt_mean only (the npz is uncompressed, so HfFileSystem seeks).
Writes data/full_<idx>.npz with W (16,1024,1024) f32, truth (16,1024), oracle mu/var/k3/k4 (16,1024).
usage: python3 fetch_full.py SHARD [MAX_MLPS]"""
import os, sys, io
import numpy as np
import pyarrow.parquet as pq
from huggingface_hub import hf_hub_download, HfFileSystem
shard = int(sys.argv[1]); mx = int(sys.argv[2]) if len(sys.argv) > 2 else 99
os.makedirs('data', exist_ok=True)
pf = hf_hub_download('aicrowd/arc-whestbench-public-2026', f'data/full-{shard:05d}-of-00067.parquet', repo_type='dataset', revision='v2-phase2')
t = pq.read_table(pf)
cols = t.column_names; print('columns', cols)
fs = HfFileSystem()
for r in range(min(mx, t.num_rows)):
    row = {c: t.column(c)[r].as_py() for c in ('mlp_id', 'mlp_name', 'mlp_seed')}
    idx = int(row['mlp_id'])
    W = np.asarray(t.column('weights')[r].as_py(), dtype=np.float32)
    truth = np.asarray(t.column('all_layer_means')[r].as_py(), dtype=np.float64)
    with fs.open(f'datasets/keenanpepper/arc-whestbench-p2-full1000-N1e9/moments/mlp_{idx:05d}.npz', 'rb') as f:
        z = np.load(f)
        m1, m2, m3, m4 = z['pre_mean'], z['pre_m2'], z['pre_m3'], z['pre_m4']
        gt = z['gt_mean']; name = str(z['mlp_name'])
    assert name == row['mlp_name'], (name, row['mlp_name'])
    var = m2 - m1 ** 2
    k3 = m3 - 3 * m1 * m2 + 2 * m1 ** 3
    k4 = m4 - 4 * m1 * m3 - 3 * m2 ** 2 + 12 * m1 ** 2 * m2 - 6 * m1 ** 4
    np.savez(f'data/full_{idx:05d}.npz', W=W, truth=truth, mu=m1, var=var, k3=k3, k4=k4, name=name, seed=row['mlp_seed'])
    print('saved', idx, name, 'rms(truth-gt)', float(np.sqrt(((truth - gt) ** 2).mean())), flush=True)
