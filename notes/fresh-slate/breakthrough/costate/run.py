"""Run co-state variants on a bench set; one JSON line per (variant, mlp) appended to results/<set>.jsonl."""
import sys, os, json, time, argparse, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "../../bench"))
import bench, costate

VARIANTS = {
    "gauss": dict(A=-1),
    "full": dict(),
    "A0": dict(A=0), "A1": dict(A=1), "A2": dict(A=2), "A3": dict(A=3), "A5": dict(A=5), "A7": dict(A=7),
    "A0slice": dict(A=0, old="slice"), "A1slice": dict(A=1, old="slice"), "A3slice": dict(A=3, old="slice"),
    "nocoinc": dict(coinc=False),
    "A0gs": dict(A=0, old="gsmslice"), "A2gs": dict(A=2, old="gsmslice"), "A1gs": dict(A=1, old="gsmslice"), "A3gs": dict(A=3, old="gsmslice"),
    "A1gl_nc": dict(A=1, old="gsm", law=True, coinc=False), "A2gl_nc": dict(A=2, old="gsm", law=True, coinc=False),
    "A3gl_nc": dict(A=3, old="gsm", law=True, coinc=False), "A3gsl_nc": dict(A=3, old="gsmslice", law=True, coinc=False),
    "A1gl": dict(A=1, old="gsm", law=True), "A2gl": dict(A=2, old="gsm", law=True), "A3gl": dict(A=3, old="gsm", law=True),
    "A1gsl": dict(A=1, old="gsmslice", law=True), "A3gsl": dict(A=3, old="gsmslice", law=True),
    "A0gp": dict(A=0, old="gsm"), "A1gp": dict(A=1, old="gsm"), "A2gp": dict(A=2, old="gsm"), "A3gp": dict(A=3, old="gsm"),
    "A1hfit": dict(A=1, oldfilter="hfit"), "A1hfitd": dict(A=1, oldfilter="hfitd"), "A0hfit": dict(A=0, oldfilter="hfit"),
    "A1oldD": dict(A=1, oldfilter="diag"), "A1gsm": dict(A=1, oldfilter="gsm"), "A1gsmoff": dict(A=1, oldfilter="gsm_off"),
    "A0gsm": dict(A=0, oldfilter="gsm"), "A0gsmoff": dict(A=0, oldfilter="gsm_off"), "A3gsm": dict(A=3, oldfilter="gsm"),
    "A0R1": dict(A=0, oldfilter="rank", r=1), "A0R8": dict(A=0, oldfilter="rank", r=8), "A1oldOff": dict(A=1, oldfilter="off"),
    "A1oldR1": dict(A=1, oldfilter="rank", r=1), "A1oldR8": dict(A=1, oldfilter="rank", r=8),
    "A1oldR64": dict(A=1, oldfilter="rank", r=64), "A1oldR256": dict(A=1, oldfilter="rank", r=256),
    "A0p0": dict(A=0, old="pool", Ap=0), "A0p1": dict(A=0, old="pool", Ap=1), "A0p3": dict(A=0, old="pool", Ap=3),
    "A1p0": dict(A=1, old="pool", Ap=0), "A1p1": dict(A=1, old="pool", Ap=1), "A1p3": dict(A=1, old="pool", Ap=3),
    "A3p0": dict(A=3, old="pool", Ap=0), "A3p3": dict(A=3, old="pool", Ap=3),
    "A0pinf": dict(A=0, old="pool", Ap=None), "A1pinf": dict(A=1, old="pool", Ap=None),
    "rk_half": dict(rank=lambda a: None if a < 1 else max(32, 1024 // (2 * a))),
}

ap = argparse.ArgumentParser()
ap.add_argument("--set", required=True); ap.add_argument("--variants", required=True)
ap.add_argument("--mlps", default=None)
a = ap.parse_args()
S = bench.load_set(a.set)
mlps = range(len(S["seeds"])) if a.mlps is None else [int(x) for x in a.mlps.split(",")]
os.makedirs(os.path.join(HERE, "results_live"), exist_ok=True)
fn = os.path.join(HERE, "results_live", f"{a.set}{os.environ.get('RES_SUFFIX', '')}.jsonl")
for i in mlps:
    W = bench.weights(S, i).astype(np.float64)
    truth = S["means"][i]
    for v in a.variants.split(","):
        kw = dict(VARIANTS[v])
        if v == "rk_half":
            n = W.shape[1]; kw["rank"] = (lambda nn: (lambda age: None if age < 1 else max(16, nn // (2 * age))))(n)
        rec = []
        t0 = time.time()
        pred = costate.predict(W, record=rec, **kw)
        dt = time.time() - t0
        lay = ((pred - truth) ** 2).mean(1)
        if os.environ.get("PRED_DIR"):
            os.makedirs(os.environ["PRED_DIR"], exist_ok=True)
            np.save(os.path.join(os.environ["PRED_DIR"], f"{a.set}_{v}_{i}.npy"), pred)
        row = dict(set=a.set, variant=v, mlp=i, mse=float(lay[-1]), raw=float(lay[-1] - S["noise"][i]),
                   noise=float(S["noise"][i]), all_layer=float(lay.mean()), per_layer=[float(x) for x in lay],
                   nprod=rec[-1]["nprod"], sec=dt)
        with open(fn, "a") as f:
            f.write(json.dumps(row) + "\n")
        print(f"{a.set} {v:9s} mlp {i}: raw {row['raw']:.3e}  nprod {row['nprod']}  {dt:.0f}s", flush=True)
