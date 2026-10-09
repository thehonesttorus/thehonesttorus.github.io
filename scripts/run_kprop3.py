"""python scripts/run_kprop3.py DATA_DIR NET [NET ...] [radial=0|1] [tag=...]  -- numpy port of the reference K=3 simple chain"""
import sys, os, time, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from whest.kprop3 import kprop3_chain
D = sys.argv[1]; nets = [int(a) for a in sys.argv[2:] if "=" not in a]; opts = dict(radial=1, tag="kprop3")
for a in sys.argv[2:]:
    if "=" in a:
        k, v = a.split("="); opts[k] = int(v) if k == "radial" else v
for net in nets:
    W = np.load(f"{D}/W_off{net}.npy"); mt = np.load(f"{D}/truth_off{net}.npz")["m"].astype(np.float64)
    t0 = time.time(); rec = {}; out = kprop3_chain(W, record=rec, radial=bool(opts["radial"])); dt = time.time() - t0
    per = np.mean((out - mt) ** 2, axis=1)
    line = f"net {net} {opts['tag']} radial={opts['radial']}: final MSE {per[-1]:.4e} ({dt:.0f}s) | per layer: " + " ".join(f"{p:.1e}" for p in per)
    print(line, flush=True)
    if os.environ.get("OUT"):
        np.savez(f"{os.environ['OUT']}/{opts['tag']}_{net}.npz", out=out, per=per, c4=np.array([rec[l]["c4"] for l in rec]),
                 K3_21=np.array([rec[l]["K3_21"] for l in rec]), K3_3=np.array([rec[l]["K3_3"] for l in rec]), var=np.array([rec[l]["var"] for l in rec]), m=np.array([rec[l]["m"] for l in rec]))
        open(f"{os.environ['OUT']}/{opts['tag']}_{net}.txt", "w").write(line + "\n")
