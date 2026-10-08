import sys, os, time, importlib.util, numpy as np, flopscope as flops
from whestbench import MLP
os.environ["K3_WIN"] = "0"
i = int(sys.argv[1])
Wcol = np.load(f"../official/W_off{i}.npy"); mt = np.load(f"../official/truth_off{i}.npz")["m"].astype(np.float64)
mlp = MLP(width=1024, depth=16, weights=[np.ascontiguousarray(W.T).astype(np.float32) for W in Wcol])
for cfg in sys.argv[2:]:
    r, age = cfg.split(":"); os.environ["K3_TRUNC_R"] = r; os.environ["K3_TRUNC_AGE"] = age
    spec = importlib.util.spec_from_file_location("est", "est_win.py"); mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    t0 = time.time()
    with flops.BudgetContext(flop_budget=2**42, wall_time_limit_s=1800.0, quiet=True) as ctx:
        out = np.asarray(mod.Estimator().predict(mlp, 2**41), dtype=np.float64)
    per = np.mean((out - mt)**2, axis=1)
    print(f"net {i} truncate old legs to rank {r} for age >= {age}: final MSE {per[-1]:.4e}  layers 8,12: {per[8]:.2e} {per[12]:.2e}  ({time.time()-t0:.0f}s)", flush=True)
