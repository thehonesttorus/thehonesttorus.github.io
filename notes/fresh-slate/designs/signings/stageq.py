"""Stage Q: run an estimator on the baked MLPs and report per-layer MSE vs truth (truth noise subtracted)."""
import sys, os, glob, json, numpy as np, importlib
TRUTH = "/root/sg/truth"

def load(width, seed):
    d = np.load(f"{TRUTH}/w{width}_s{seed}.npz"); W = np.load(f"{TRUTH}/W_w{width}_s{seed}.npy").astype(float)
    N = int(d["n"]); mean = d["S"] / N; var = d["Q"] / N - mean ** 2
    return W, mean, var / N, N

def evaluate(est_fn, widths, seeds, tag=""):
    rows = []
    for w in widths:
        for s in seeds:
            if not os.path.exists(f"{TRUTH}/w{w}_s{s}.npz"): continue
            W, mean, noise, N = load(w, s)
            est = est_fn(W)
            err = ((est - mean) ** 2).mean(1) - noise.mean(1)
            rows.append(dict(width=w, seed=s, N=N, final=float(err[-1]), noise=float(noise[-1].mean()),
                             layers=[float(e) for e in err]))
            print(tag, w, s, f"N={N:.2e} final raw MSE {err[-1]:.3e} (noise {noise[-1].mean():.1e})  L2 {err[1]:.2e} L4 {err[3]:.2e} L8 {err[7]:.2e}", flush=True)
    return rows

if __name__ == "__main__":
    mod, fn = sys.argv[1].split(":"); cfg = json.loads(sys.argv[2]) if len(sys.argv) > 2 else {}
    widths = [int(x) for x in (sys.argv[3] if len(sys.argv) > 3 else "64,128,256").split(",")]
    est = getattr(importlib.import_module(mod), fn)
    rows = evaluate(lambda W: est(W, cfg), widths, [int(x) for x in (sys.argv[5] if len(sys.argv) > 5 else "0,1,2,3").split(",")], tag=sys.argv[1] + json.dumps(cfg))
    out = sys.argv[4] if len(sys.argv) > 4 else None
    if out: json.dump(rows, open(out, "w"), indent=1)
