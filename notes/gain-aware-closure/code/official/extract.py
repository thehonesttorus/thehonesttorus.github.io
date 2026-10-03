# Extract official WhestBench Phase-2 MLPs: weights in column convention (h_l = relu(W h_{l-1}), i.e. W = W_official^T)
import pyarrow.parquet as pq, numpy as np, sys
f = pq.ParquetFile(sys.argv[1])
for rg in range(f.num_row_groups):
    t = f.read_row_group(rg); R = t.num_rows
    W = t.column("weights").combine_chunks().flatten().flatten().flatten().to_numpy().reshape(R, 16, 1024, 1024)
    M = t.column("all_layer_means").combine_chunks().flatten().flatten().to_numpy().reshape(R, 16, 1024)
    for r in range(R):
        mid = t.column("mlp_id")[r].as_py()
        np.save(f"W_off{mid}.npy", np.ascontiguousarray(W[r].transpose(0, 2, 1)))
        np.savez(f"truth_off{mid}.npz", m=M[r], avg_variance=t.column("avg_variance")[r].as_py(), seed=t.column("mlp_seed")[r].as_py())
        print(mid, end=" ", flush=True)
