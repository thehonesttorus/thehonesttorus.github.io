"""python scripts/run_gauss.py DATA_DIR NET [NET ...] [--K k]: Gaussian closure vs truth, per layer."""
import sys, os, time, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from whest.relu_gauss import gauss_chain
D = sys.argv[1]; K = 10
args = [a for a in sys.argv[2:]]
if "--K" in args:
    K = int(args[args.index("--K") + 1]); del args[args.index("--K"):args.index("--K") + 2]
for net in map(int, args):
    W = np.load(f"{D}/W_off{net}.npy"); mt = np.load(f"{D}/truth_off{net}.npz")["m"].astype(np.float64)
    t0 = time.time(); out = gauss_chain(W, K=K); dt = time.time() - t0
    per = np.mean((out - mt) ** 2, axis=1)
    print(f"net {net} K={K} final MSE {per[-1]:.4e} ({dt:.1f}s) | per layer: " + " ".join(f"{p:.1e}" for p in per), flush=True)
