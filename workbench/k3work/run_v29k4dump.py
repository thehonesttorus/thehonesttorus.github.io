# Production chain (current best flags) with the layer dump on: keeps the post-activation (2,2) slice K22, the
# post-activation kappa4 diagonal K4v, the variance K2v, and the pre-activation g4row / mu / var at the dump layers.
import sys, os, importlib.util, numpy as np, pickle, flopscope as flops
from whestbench import MLP
i = int(sys.argv[1]); layers = sys.argv[2]
os.environ.update({"V29_WARM_JOIN": "1", "V29_WARM_FB": "1", "V17_R_RES": "4", "V26_STRASSEN": "6", "V26_STRASSEN_MIN": "16",
                   "V32_JOIN_SMM": "1", "V32_ROT_SMM": "1", "V32_JOIN_POST": "3", "V32_JP_C": "0.1", "V21_R_OLD": "320",
                   "V24_R_OLD2": "192", "V18_R_FB": "2", "V24_AGE_OLD2": "8", "V33_K4Q": "3", "V33_K4Q_RANK": "4", "V29_DUMP_LAYERS": layers})
spec = importlib.util.spec_from_file_location("estv29", "est_v29.py"); mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
Wcol = np.load(f"../official/W_off{i}.npy"); mt = np.load(f"../official/truth_off{i}.npz")["m"].astype(np.float64)
mlp = MLP(width=1024, depth=16, weights=[np.ascontiguousarray(W.T).astype(np.float32) for W in Wcol])
with flops.BudgetContext(flop_budget=2**44, wall_time_limit_s=900.0, quiet=True):
    out = np.asarray(mod.Estimator().predict(mlp, 2**41), dtype=np.float64)
keep = ("K22", "K4v", "K2v", "g4row", "mu", "var", "D21", "D3", "C_off", "wk4m")
dumps = [dict(layer=d["layer"], **{k: d[k] for k in keep if k in d}) for d in mod.DUMPS]
print("raw", np.mean((out[-1] - mt[-1])**2), "dumps", [d["layer"] for d in dumps], flush=True)
pickle.dump(dict(dumps=dumps, out=out, truth=mt), open(f"v29k4_off{i}.pkl", "wb"), protocol=4)
