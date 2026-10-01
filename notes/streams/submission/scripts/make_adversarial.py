# Adversarial robustness sets: copy a baked 1-MLP dataset, scale its weights (x10) or take |W|,
# and replace the truth with a crude float64 MC estimate (2000 samples). usage: SRC DST {x10,abs}
import sys, os, shutil, json, glob
import numpy as np, pyarrow.parquet as pq, pyarrow as pa
src, dst, mode = sys.argv[1:4]
shutil.copytree(src, dst)
f = glob.glob(os.path.join(dst, "data", "*.parquet"))[0]
t = pq.read_table(f)
sch = t.schema
w = np.array(t.column("weights").to_pylist()[0], dtype=np.float32)
if mode == "x10": w2 = w * 10
elif mode == "abs": w2 = np.abs(w)
elif mode == "x0.1": w2 = w * 0.1
else: raise SystemExit(mode)
# crude MC truth so the MSE column is meaningful-ish (float64, 2000 samples)
rng = np.random.default_rng(0)
x = rng.standard_normal((2000, w.shape[1])).astype(np.float64)
means = []
for l in range(w2.shape[0]):
    x = np.maximum(x @ w2[l].astype(np.float64), 0.0); means.append(x.mean(0))
means = np.array(means, dtype=np.float32)
print(mode, "final mean magnitude", float(np.abs(means[-1]).mean()), "finite", bool(np.isfinite(means).all()))
cols = {n: t.column(n) for n in t.column_names}
cols["weights"] = pa.array([w2.tolist()], type=sch.field("weights").type)
cols["all_layer_means"] = pa.array([means.tolist()], type=sch.field("all_layer_means").type)
cols["final_means"] = pa.array([means[-1].tolist()], type=sch.field("final_means").type)
pq.write_table(pa.table(cols, schema=sch), f)
