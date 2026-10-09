"""python scripts/run_k3v2.py DATA_DIR NET [NET ...] [opt=val ...]"""
import sys, os, time, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from whest.k3chain3 import k3_chain3 as k3_chain2
D = sys.argv[1]; nets = [int(a) for a in sys.argv[2:] if "=" not in a]; opts = {}
for a in sys.argv[2:]:
    if "=" in a:
        k, v = a.split("="); opts[k] = v if k in ("k4", "k22gate") else (float(v) if "." in v else int(v))
tag = os.environ.get("TAG", "k3v2")
for net in nets:
    W = np.load(f"{D}/W_off{net}.npy"); mt = np.load(f"{D}/truth_off{net}.npz")["m"].astype(np.float64)
    t0 = time.time(); rec = {}; out, fl = k3_chain2(W, opts, record=rec); dt = time.time() - t0
    per = np.mean((out - mt) ** 2, axis=1)
    line = f"net {net} {tag} {opts}: final MSE {per[-1]:.4e}  est C/B {fl / 2 ** 41:.3f}  ({dt:.0f}s) | per layer: " + " ".join(f"{p:.1e}" for p in per)
    print(line, flush=True)
    if os.environ.get("OUT"):
        np.savez(f"{os.environ['OUT']}/{tag}_{net}.npz", out=out, per=per, flops=fl, D3=np.array([rec[l]["D3"] for l in rec]))
        open(f"{os.environ['OUT']}/{tag}_{net}.txt", "w").write(line + "\n")
