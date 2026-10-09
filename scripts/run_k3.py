"""python scripts/run_k3.py DATA_DIR NET [NET ...] [opt=val ...]   e.g. fold=1 hub=1 k4=path window=0"""
import sys, os, time, json, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from whest.k3chain import k3_chain
D = sys.argv[1]; nets = [int(a) for a in sys.argv[2:] if "=" not in a]
opts = {}
for a in sys.argv[2:]:
    if "=" in a:
        k, v = a.split("="); opts[k] = v if k in ("k4",) else (float(v) if "." in v else int(v))
tag = os.environ.get("TAG", "k3")
for net in nets:
    W = np.load(f"{D}/W_off{net}.npy"); mt = np.load(f"{D}/truth_off{net}.npz")["m"].astype(np.float64)
    t0 = time.time(); rec = {}; out, fl = k3_chain(W, opts, record=rec); dt = time.time() - t0
    per = np.mean((out - mt) ** 2, axis=1)
    line = f"net {net} {tag} {opts}: final MSE {per[-1]:.4e}  est C/B {fl / 2 ** 41:.3f}  ({dt:.0f}s) | per layer: " + " ".join(f"{p:.1e}" for p in per)
    print(line, flush=True)
    if os.environ.get("OUT"):
        np.savez(f"{os.environ['OUT']}/{tag}_{net}.npz", out=out, per=per, flops=fl, D3=np.array([rec[l]["D3"] for l in rec]), K4=np.array([rec[l]["K4"] for l in rec]))
        open(f"{os.environ['OUT']}/{tag}_{net}.txt", "w").write(line + "\n")
