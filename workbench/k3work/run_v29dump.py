import sys, os, importlib.util, numpy as np, pickle, flopscope as flops
from whestbench import MLP
i = int(sys.argv[1]); layers = sys.argv[2]
os.environ["V29_WARM_JOIN"] = "1"; os.environ["V29_WARM_FB"] = "1"; os.environ["V29_DUMP_LAYERS"] = layers
spec = importlib.util.spec_from_file_location("estv29", "est_v29.py"); mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
Wcol = np.load(f"../official/W_off{i}.npy"); mt = np.load(f"../official/truth_off{i}.npz")["m"].astype(np.float64)
mlp = MLP(width=1024, depth=16, weights=[np.ascontiguousarray(W.T).astype(np.float32) for W in Wcol])
with flops.BudgetContext(flop_budget=2**41, wall_time_limit_s=900.0, quiet=True):
    out = np.asarray(mod.Estimator().predict(mlp, 2**41), dtype=np.float64)
print("raw", np.mean((out[-1]-mt[-1])**2), "dumps", len(mod.DUMPS), flush=True)
pickle.dump(mod.DUMPS, open(f"v29dump_off{i}.pkl", "wb"))
