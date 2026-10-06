# K4 diagonal tests: K4_DERIVED=1 (3 g var^2 from the chain's own D3), =2 (log only), or G4SCALE on all layers.
import sys, os, importlib.util, numpy as np, pickle, flopscope as flops
from whestbench import MLP
i = int(sys.argv[1]); mode = sys.argv[2]; scale = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
os.environ["K3_WIN"] = "0"
if mode == "derived": os.environ["K4_DERIVED"] = "1"; os.environ["K4_DERIVED_SCALE"] = str(scale)
elif mode == "log": os.environ["K4_DERIVED"] = "2"
spec = importlib.util.spec_from_file_location("est", "est_win.py"); mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
if mode == "scale": mod.G4SCALE = scale; mod.G4LAYERS = tuple(range(16))
Wcol = np.load(f"../official/W_off{i}.npy"); mt = np.load(f"../official/truth_off{i}.npz")["m"].astype(np.float64)
mlp = MLP(width=1024, depth=16, weights=[np.ascontiguousarray(W.T).astype(np.float32) for W in Wcol])
with flops.BudgetContext(flop_budget=2**42, wall_time_limit_s=900.0, quiet=True):
    out = np.asarray(mod.Estimator().predict(mlp, 2**41), dtype=np.float64)
print(f"net {i} mode={mode} scale={scale}: final MSE {np.mean((out[-1]-mt[-1])**2):.4e}", flush=True)
if mod.K4LOG:
    for r in mod.K4LOG: print(f"   layer {r['layer']:2d}: g {r['g']:.5f} (corr {r['corr']:+.3f})  chain g4row = {r['old_over_var2']:.5f} var^2 (corr {r['old_corr_var2']:+.3f})  3g = {3*r['g']:.5f}  ratio chain/derived {r['old_over_var2']/(3*r['g']):.3f}")
