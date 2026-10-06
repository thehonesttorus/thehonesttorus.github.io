import sys, os, importlib.util, numpy as np, pickle, flopscope as flops
from whestbench import MLP
os.environ["K3_WIN"] = "0"; os.environ["K3_SPEC"] = "1"
spec = importlib.util.spec_from_file_location("est", "est_win.py"); mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
i = int(sys.argv[1])
Wcol = np.load(f"../official/W_off{i}.npy")
mlp = MLP(width=1024, depth=16, weights=[np.ascontiguousarray(W.T).astype(np.float32) for W in Wcol])
with flops.BudgetContext(flop_budget=2**42, wall_time_limit_s=900.0, quiet=True):
    out = np.asarray(mod.Estimator().predict(mlp, 2**41), dtype=np.float64)
pickle.dump(mod.SPEC, open(f"spec_off{i}.pkl", "wb"))
for d in mod.SPEC:
    sa, sp = d["sA"], d["sP"]
    ea = (sa.sum()**2) / (sa**2).sum(); ep = (sp.sum()**2) / (sp**2).sum()
    fa = lambda r: (sa[:r]**2).sum() / (sa**2).sum(); fp = lambda r: (sp[:r]**2).sum() / (sp**2).sum()
    print(f"layer {d['layer']:2d} source born {d['source']:2d} age {d['age']:2d}: A s1/s64/s256 {sa[0]:.2e}/{sa[63]:.2e}/{sa[255]:.2e} energy@64 {fa(64):.3f} @256 {fa(256):.3f} | P s1/s64/s256 {sp[0]:.2e}/{sp[63]:.2e}/{sp[255]:.2e} energy@64 {fp(64):.3f} @256 {fp(256):.3f}", flush=True)
