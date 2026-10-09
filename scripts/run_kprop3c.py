"""python scripts/run_kprop3c.py DATA_DIR NET [window=2] [k=128] [c4scale=1.0] [tag=...]   -- two-tier K=3 chain"""
import sys, os, time, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from whest.kprop3c import kprop3c_chain
D = sys.argv[1]; nets = [int(a) for a in sys.argv[2:] if "=" not in a]; opts = {}; tag = "kprop3c"; save = 0
for a in sys.argv[2:]:
    if "=" in a:
        k, v = a.split("=")
        if k == "tag": tag = v
        elif k == "save": save = int(v)
        elif k in ("oldmode", "tier", "specres"): opts[k] = v
        else: opts[k] = float(v) if "." in v else int(v)
for net in nets:
    W = np.load(f"{D}/W_off{net}.npy"); mt = np.load(f"{D}/truth_off{net}.npz")["m"].astype(np.float64)
    t0 = time.time(); out = kprop3c_chain(W, opts); dt = time.time() - t0; per = np.mean((out - mt) ** 2, axis=1)
    line = f"net {net} {tag} {opts}: final MSE {per[-1]:.4e} ({dt:.0f}s) | per layer: " + " ".join(f"{p:.1e}" for p in per)
    print(line, flush=True)
    if os.environ.get("OUT"): open(f"{os.environ['OUT']}/{tag}_{net}.txt", "a").write(line + "\n")
    if save and os.environ.get("OUT"):
        np.save(f"{os.environ['OUT']}/{tag}_{net}_" + "_".join(f"{k}{v}" for k, v in sorted(opts.items())) + ".npy", out[-1])
