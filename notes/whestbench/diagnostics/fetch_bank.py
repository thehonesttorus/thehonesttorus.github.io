# Held-out bank: seeds + 1e9-sample ground truth for all 100 networks of the mini split, by HTTP range
# reads of the small columns only.  Weights are regenerated from mlp_seed (verified against the exact
# layer-1 means in notes/checks/phase2_output.txt).
import fsspec, numpy as np, pyarrow as pa, pyarrow.parquet as pq, time
BASE=("https://huggingface.co/datasets/aicrowd/arc-whestbench-public-2026/resolve/v2-phase2/data/"
      "mini-{:05d}-of-00007.parquet")
def to_np(col):
    a=col.combine_chunks() if isinstance(col,pa.ChunkedArray) else col
    while pa.types.is_list(a.type) or pa.types.is_fixed_size_list(a.type) or pa.types.is_large_list(a.type): a=a.flatten()
    return a.to_numpy(zero_copy_only=False)
fs=fsspec.filesystem("https"); S=[]; Y=[]; N=[]; t0=time.time()
for i in range(7):
    for attempt in range(4):
        try:
            with fs.open(BASE.format(i), block_size=2**20) as f:
                tb=pq.ParquetFile(f).read_row_group(0, columns=["mlp_name","mlp_seed","all_layer_means"])
            break
        except Exception as e:
            print(f"  shard {i} attempt {attempt}: {e}", flush=True); time.sleep(2**attempt)
    else: raise SystemExit(f"shard {i} failed")
    S += tb.column("mlp_seed").to_pylist(); N += tb.column("mlp_name").to_pylist()
    Y.append(to_np(tb.column("all_layer_means")).astype(np.float32).reshape(tb.num_rows,16,1024))
    print(f"shard {i}: {tb.num_rows} nets, {time.time()-t0:.0f}s", flush=True)
Y=np.concatenate(Y,0)
np.savez_compressed("bank.npz", seeds=np.array(S,dtype=np.int64), Y=Y, names=np.array(N))
print("bank:", Y.shape, f"{time.time()-t0:.0f}s")
