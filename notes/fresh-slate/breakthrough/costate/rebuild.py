"""Rebuild result rows from saved predictions (PRED_DIR) against bench truth; nprod/sec taken from the run log."""
import sys, os, json, re, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "../../bench"))
import bench
P = os.environ["PRED_DIR"]; name, variants, log, out = sys.argv[1], sys.argv[2].split(","), sys.argv[3], sys.argv[4]
S = bench.load_set(name)
meta = {}
for line in open(log):
    m = re.match(r"\S+ (\S+)\s+mlp (\d+): raw \S+\s+nprod (\d+)\s+(\d+)s", line)
    if m: meta[(m.group(1), int(m.group(2)))] = (int(m.group(3)), float(m.group(4)))
rows = []
for v in variants:
    for i in range(len(S["seeds"])):
        f = os.path.join(P, f"{name}_{v}_{i}.npy")
        if not os.path.exists(f) or (v, i) not in meta: continue
        pred = np.load(f); lay = ((pred - S["means"][i]) ** 2).mean(1)
        rows.append(dict(set=name, variant=v, mlp=i, mse=float(lay[-1]), raw=float(lay[-1] - S["noise"][i]), noise=float(S["noise"][i]),
                         all_layer=float(lay.mean()), per_layer=[float(x) for x in lay], nprod=meta[(v, i)][0], sec=meta[(v, i)][1]))
with open(out, "w") as fo:
    for r in rows: fo.write(json.dumps(r) + "\n")
print(len(rows), "rows")
