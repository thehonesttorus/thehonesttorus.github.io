# K3 chain (shape-generic base path) on the small-width nets, for reference against the arrow mixture.
import sys, os, importlib.util, numpy as np, flopscope as flops
from whestbench import MLP
os.environ["K3_WIN"] = "0"
spec = importlib.util.spec_from_file_location("est", "est_win.py"); mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
n, L = int(sys.argv[1]), int(sys.argv[2])
for tag in sys.argv[3:]:
    Wc = np.load(f"../num13/W_n{n}_L{L}_s{tag}.npy") if n == 64 else np.load(f"../num12/W_n{n}_L{L}_s{tag}.npy")
    mt = (np.load(f"../num13/truth_n{n}_L{L}_s{tag}.npz") if n == 64 else np.load(f"../num12/truth_n{n}_L{L}_s{tag}.npz"))["m"]
    mlp = MLP(width=n, depth=L, weights=[np.ascontiguousarray(W.T).astype(np.float32) for W in Wc])
    with flops.BudgetContext(flop_budget=2**42, wall_time_limit_s=900.0, quiet=True) as ctx:
        out = np.asarray(mod.Estimator().predict(mlp, 2**41), dtype=np.float64)
    print(f"n={n} L={L} s={tag}: K3 chain final MSE {np.mean((out[-1]-mt[-1])**2):.3e}  flops {ctx.flops_used/2**41:.4f} B", flush=True)
