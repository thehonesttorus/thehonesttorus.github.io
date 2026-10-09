# The adopted system's own two-site slices per layer (note XLIII section 3): what its Wick stage reads after the V47
# calibration. python chaindump2.py NET  ->  $OUT/chaindump2_{NET}.npz with, per layer l: D3_l, D21_l, g4row_l, wk4m_l,
# wk431_l (wk431[a, c] = kappa(z_a, z_c, z_c, z_c)), var_l, C_off_l, and out (16 x n predicted means).
import os, sys, importlib.util, numpy as np, flopscope as flops
from whestbench import MLP
net = int(sys.argv[1])
os.environ["V29_DUMP_LAYERS"] = ",".join(str(x) for x in range(16))
spec = importlib.util.spec_from_file_location("estf", "estimator_final_v56.py"); mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
OFF = os.environ.get("OFFDIR", "../official")
Wcol = np.load(f"{OFF}/W_off{net}.npy")
mlp = MLP(width=1024, depth=16, weights=[np.ascontiguousarray(W.T).astype(np.float32) for W in Wcol])
with flops.BudgetContext(flop_budget=2**45, wall_time_limit_s=3000.0, quiet=True):
    out = np.asarray(mod.Estimator().predict(mlp, 2**41), dtype=np.float64)
mt = np.load(f"{OFF}/truth_off{net}.npz")["m"].astype(np.float64)
keep = ("D3", "D21", "g4row", "wk4m", "wk431", "var", "C_off")
got = {}
for d in mod.DUMPS:
    if "wk4m" not in d:
        continue
    for k in keep:
        if d.get(k) is not None and f"{k}_{d['layer']}" not in got:
            got[f"{k}_{d['layer']}"] = np.asarray(d[k], np.float32)
print(f"net {net}: final raw {np.mean((out[-1] - mt[-1]) ** 2):.4e}; layers with slices "
      f"{sorted(int(k.split('_')[-1]) for k in got if k.startswith('D21_'))}", flush=True)
np.savez(f"{os.environ.get('OUT', '.')}/chaindump2_{net}.npz", out=out, **got)
