"""Dilation-sector experiment inside FC (region stream) on w1024_d16. Appends rows to results_live/fcdil_w1024.jsonl and
the full-run instrumentation to results_live/fcdil_instr_<mlp>.json."""
import sys, os, json, time, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "../../bench")); sys.path.insert(0, os.path.join(HERE, "../region"))
import bench, fc_dil
S = bench.load_set("w1024_d16")
mlps = [int(x) for x in sys.argv[1].split(",")]
variants = sys.argv[2].split(",") if len(sys.argv) > 2 else ["full", "drop2", "drop4", "oracle2", "oracle4", "gl2", "gl4", "gp2", "gp4"]
os.makedirs(os.path.join(HERE, "results_live"), exist_ok=True)
for i in mlps:
    W = bench.weights(S, i); T = S["means"][i]
    for v in variants:
        kw = dict(slices=2, k4mf=True)
        inst = None
        if v == "full":
            inst = []; kw["instr"] = inst
        else:
            mode, A = v[:-1], int(v[-1])
            kw.update(dil=mode, A=A)
        t0 = time.time(); p = fc_dil.run(W, **kw); dt = time.time() - t0
        raw = float(((p[-1] - T[-1]) ** 2).mean() - S["noise"][i])
        row = dict(mlp=i, variant=v, raw=raw, sec=dt, layers=((p - T) ** 2).mean(1).tolist())
        open(os.path.join(HERE, "results_live", "fcdil_w1024.jsonl"), "a").write(json.dumps(row) + "\n")
        if inst is not None:
            json.dump(inst, open(os.path.join(HERE, "results_live", f"fcdil_instr_{i}.json"), "w"))
        print(f"mlp {i} {v:8s} raw {raw:.3e}  {dt:.0f}s", flush=True)
