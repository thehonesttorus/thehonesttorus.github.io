# K3 chain wrapped in GAC (gain-conditional state, E[G] output) on official networks; env: GAIN_COND, GAIN_NO_G4, V17_LAM_SCALE, GAIN_SCALE
import sys, os, time, importlib.util, numpy as np, flopscope as flops
from whestbench import MLP
os.environ.setdefault("K3_WIN", "0")
spec = importlib.util.spec_from_file_location("est", "est_win.py"); mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
tag = f"GAIN_COND={os.environ.get('GAIN_COND','0')} NO_G4={os.environ.get('GAIN_NO_G4','0')} LAM_SCALE={os.environ.get('V17_LAM_SCALE','1')} GAIN_SCALE={os.environ.get('GAIN_SCALE','1')} WIN={os.environ.get('K3_WIN','0')}"
for i in [int(x) for x in sys.argv[1].split(",")]:
    Wcol = np.load(f"../official/W_off{i}.npy"); mt = np.load(f"../official/truth_off{i}.npz")["m"].astype(np.float64)
    mlp = MLP(width=1024, depth=16, weights=[np.ascontiguousarray(W.T).astype(np.float32) for W in Wcol])
    t0 = time.time()
    with flops.BudgetContext(flop_budget=2**42, wall_time_limit_s=900.0, quiet=True) as ctx:
        out = np.asarray(mod.Estimator().predict(mlp, 2**41), dtype=np.float64)
    per = np.mean((out - mt)**2, axis=1)
    print(f"net {i} [{tag}]: final MSE {per[-1]:.4e}  C/B {ctx.flops_used/2**41:.4f}  {time.time()-t0:.0f}s  layers 4,8,12,15: {per[4]:.2e} {per[8]:.2e} {per[12]:.2e} {per[15]:.2e}", flush=True)
