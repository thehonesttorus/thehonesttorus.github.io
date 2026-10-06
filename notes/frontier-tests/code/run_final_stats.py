# Compare the chain's final-layer readout inputs (var, D3, kappa4 diagonal) with 524k-sample Monte Carlo moments of z16.
import sys, os, importlib.util, numpy as np, flopscope as flops
from whestbench import MLP
from scipy.special import ndtr
os.environ["K3_WIN"] = "0"; os.environ["K3_DUMP2"] = "1"
spec = importlib.util.spec_from_file_location("est", "est_win.py"); mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
Wcol = np.load("../official/W_off0.npy"); mt = np.load("../official/truth_off0.npz")["m"].astype(np.float64)
mod.ORACLE_PREV = mt.astype(np.float32); mod.ORACLE_LAYERS = (15,)      # true mean entering the last layer
mlp = MLP(width=1024, depth=16, weights=[np.ascontiguousarray(W.T).astype(np.float32) for W in Wcol])
with flops.BudgetContext(flop_budget=2**42, wall_time_limit_s=900.0, quiet=True):
    out = np.asarray(mod.Estimator().predict(mlp, 2**41), dtype=np.float64)
d = mod.DUMP2[-1]
z = np.load("../readout/z16_moments_off0.npz"); S, N = z["S"], float(z["N"])
m1, m2, m3, m4 = S[1:]/N; var = m2 - m1**2; k3 = m3 - 3*m1*m2 + 2*m1**3; k4 = m4 - 4*m1*m3 + 6*m1**2*m2 - 3*m1**4 - 3*var**2
mu_true = Wcol[15].astype(np.float64) @ mt[14]
sig = np.sqrt(var); a = mu_true/sig
print(f"final MSE with true incoming mean: {np.mean((out[-1]-mt[-1])**2):.3e}")
print(f"var:  chain/MC mean ratio {np.mean(d['var']/var):.5f} (MC noise of this average ~{np.sqrt(2/N)/32:.1e}); per-neuron rel rms {np.sqrt(np.mean((d['var']/var-1)**2)):.4f} (MC noise {np.sqrt(2/N):.4f})")
print(f"D3:   regression slope chain~MC {np.sum(d['D3']*k3)/np.sum(k3*k3):.4f}; mean (chain-MC)/sig^3 {np.mean((d['D3']-k3)/sig**3):+.5f} (noise ~{np.sqrt(15/N)/32:.1e}); rms diff/sig^3 {np.sqrt(np.mean(((d['D3']-k3)/sig**3)**2)):.4f} (MC noise ~{np.sqrt(15/N):.4f})")
if d['g4row'] is not None:
    for sc in [1.0, 0.5, 2.0, 1/24, 24.0]:
        g = d['g4row']*sc
        print(f"k4 diag (g4row x {sc:g}): slope chain~MC {np.sum(g*k4)/np.sum(k4*k4):.4f}; mean (chain-MC)/sig^4 {np.mean((g-k4)/sig**4):+.5f} (noise ~{np.sqrt(96/N)/32:.1e}); rms {np.sqrt(np.mean(((g-k4)/sig**4)**2)):.4f} (MC noise ~{np.sqrt(96/N):.4f})")
np.savez("final_stats_off0.npz", var_c=d['var'], D3_c=d['D3'], g4_c=d['g4row'] if d['g4row'] is not None else 0, var_mc=var, k3_mc=k3, k4_mc=k4, mu_true=mu_true)
