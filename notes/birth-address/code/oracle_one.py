# One chain run for the oracle attribution (note XXXI).   python oracle_one.py NET ORACLE MCFILE TAG
# ORACLE: none | dump | D3 | D21 | D3+D21. "dump" saves the chain's own per-layer D3, D21, var for comparison.
import sys, os, importlib.util, numpy as np, flopscope as flops
from whestbench import MLP
net, orc, mcf, tag = int(sys.argv[1]), sys.argv[2], sys.argv[3], sys.argv[4]
os.environ.update({"V29_WARM_JOIN": "1", "V29_WARM_FB": "1", "V17_R_RES": "4", "V26_STRASSEN": "6", "V26_STRASSEN_MIN": "16",
                   "V32_JOIN_SMM": "1", "V32_ROT_SMM": "1", "V32_JOIN_POST": "3", "V32_JP_C": "0.1", "V21_R_OLD": "320",
                   "V24_R_OLD2": "192", "V18_R_FB": "2", "V24_AGE_OLD2": "8", "V33_K4Q": "3", "V33_K4Q_RANK": "4",
                   "V34_OPT": "abcde"})
if orc == "dump":
    os.environ["V29_DUMP_LAYERS"] = ",".join(str(x) for x in range(16))
elif orc != "none":
    os.environ["V37_ORACLE"] = orc.replace("+", ","); os.environ["V37_ORACLE_FILE"] = mcf
spec = importlib.util.spec_from_file_location("estv29", "est_v29.py"); mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
Wcol = np.load(f"../official/W_off{net}.npy"); mt = np.load(f"../official/truth_off{net}.npz")["m"].astype(np.float64)
mlp = MLP(width=1024, depth=16, weights=[np.ascontiguousarray(W.T).astype(np.float32) for W in Wcol])
with flops.BudgetContext(flop_budget=2**44, wall_time_limit_s=900.0, quiet=True):
    out = np.asarray(mod.Estimator().predict(mlp, 2**41), dtype=np.float64)
print(f"net {net} {tag}: raw {np.mean((out[-1] - mt[-1]) ** 2):.5e}", flush=True)
if orc == "dump":
    d3 = {}; d21 = {}; var = {}
    for d in mod.DUMPS:
        if d.get("D3") is not None and d["layer"] not in d3:
            d3[d["layer"]] = np.asarray(d["D3"], np.float32); var[d["layer"]] = np.asarray(d["var"], np.float32)
        if d.get("D21") is not None and d["layer"] not in d21:
            d21[d["layer"]] = np.asarray(d["D21"], np.float32)
    np.savez(f"chain_off{net}.npz", **{f"D3_{k}": v for k, v in d3.items()}, **{f"D21_{k}": v for k, v in d21.items()},
             **{f"var_{k}": v for k, v in var.items()})
