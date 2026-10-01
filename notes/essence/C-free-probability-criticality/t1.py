"""T1: leak-corrected dilation charge. python t1.py <set> <variants> [mlps]
variant name: <base>@<leakfac>, base in costate run.py VARIANTS (e.g. A3gsl_nc@1.7). Rows -> results/t1_<set>.jsonl"""
import sys, os, json, time, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "../../fresh-slate/bench"))
import bench, costate_c
BASE = {"A1gl_nc": dict(A=1, old="gsm", law=True, coinc=False), "A2gl_nc": dict(A=2, old="gsm", law=True, coinc=False),
        "A3gl_nc": dict(A=3, old="gsm", law=True, coinc=False), "A3gsl_nc": dict(A=3, old="gsmslice", law=True, coinc=False),
        "A1gsl_nc": dict(A=1, old="gsmslice", law=True, coinc=False), "full_nc": dict(coinc=False), "gauss": dict(A=-1), "Tgsl_nc": dict(A=None, old="gsmslice", law=True, coinc=False), "Tgl_nc": dict(A=None, old="gsm", law=True, coinc=False)}
name, variants = sys.argv[1], sys.argv[2].split(",")
S = bench.load_set(name)
mlps = range(len(S["seeds"])) if len(sys.argv) < 4 else [int(x) for x in sys.argv[3].split(",")]
fn = os.path.join(HERE, "results", f"t1_{name}.jsonl")
for i in mlps:
    W = bench.weights(S, i).astype(np.float64); truth = S["means"][i]
    for v in variants:
        b, lf = (v.split("@") + ["1"])[:2]
        tau = None
        if "~" in b: b, tau = b.split("~"); tau = float(tau)
        rec = []; t0 = time.time()
        pred = costate_c.predict(W, record=rec, leakfac=float(lf), tau=tau, **BASE[b])
        lay = ((pred - truth) ** 2).mean(1)
        row = dict(set=name, variant=v, mlp=i, raw=float(lay[-1] - S["noise"][i]), nprod=rec[-1]["nprod"], sec=time.time() - t0)
        open(fn, "a").write(json.dumps(row) + "\n")
        print(f"{name} {v:14s} mlp {i}: raw {row['raw']:.3e} nprod {row['nprod']} {row['sec']:.0f}s", flush=True)
