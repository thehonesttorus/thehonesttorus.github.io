# Antithetic address-modulated gating: run the chain with gate Phi + phi*g and Phi - phi*g (same g), average.
import sys, os, subprocess, numpy as np
net, seed = int(sys.argv[1]), sys.argv[2]
mt = np.load(f"../official/truth_off{net}.npz")["m"].astype(np.float64)
outs = {}
for rep in ["0", "1", "-1"]:
    f = f"pred_rep_off{net}_s{seed}_r{rep}.npy"
    if not os.path.exists(f):
        code = ("import os,numpy as np,importlib.util,flopscope as flops\nfrom whestbench import MLP\n"
                "spec=importlib.util.spec_from_file_location('est','est_win.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)\n"
                f"Wcol=np.load('../official/W_off{net}.npy')\n"
                "mlp=MLP(width=1024,depth=16,weights=[np.ascontiguousarray(W.T).astype(np.float32) for W in Wcol])\n"
                "with flops.BudgetContext(flop_budget=2**42,wall_time_limit_s=900.0,quiet=True):\n"
                f"    np.save('{f}',np.asarray(mod.Estimator().predict(mlp,2**41),dtype=np.float64))\n")
        env = dict(os.environ, K3_WIN="0", REPLICA=rep, REPLICA_SEED=seed, OMP_NUM_THREADS="1")
        subprocess.run([sys.executable, "-c", code], env=env, check=True, stderr=subprocess.DEVNULL)
    outs[rep] = np.load(f)
mse = lambda y: np.mean((y[-1] - mt[-1])**2)
avg = 0.5 * (outs["1"] + outs["-1"])
print(f"net {net} seed {seed}: plain {mse(outs['0']):.4e} | +g {mse(outs['1']):.4e}  -g {mse(outs['-1']):.4e} | antithetic mean {mse(avg):.4e}"
      f" | corr(avg-plain, truth-plain) {np.corrcoef(avg[-1]-outs['0'][-1], mt[-1]-outs['0'][-1])[0,1]:+.3f}", flush=True)
