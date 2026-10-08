import sys, os, importlib.util, numpy as np, flopscope as flops
from whestbench import MLP
os.environ.update({k: v for k, v in (kv.split("=") for kv in sys.argv[2].split())})
os.environ["V33_SAT_LOG"] = "1"
spec = importlib.util.spec_from_file_location("e", "est_v29.py"); mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
i = int(sys.argv[1]); Wcol = np.load(f"../official/W_off{i}.npy")
mlp = MLP(width=1024, depth=16, weights=[np.ascontiguousarray(W.T).astype(np.float32) for W in Wcol])
with flops.BudgetContext(flop_budget=2**44, quiet=True):
    mod.Estimator().predict(mlp, 2**41)
print("SAT", os.environ["V33_SAT"], "dropped fraction by layer:", " ".join(f"{l}:{f:.2f}" for l, f in mod.SAT_LOG))
