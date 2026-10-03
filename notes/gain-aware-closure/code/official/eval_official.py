# Official WhestBench Phase-2 networks (1e9-sample truth): closure, closure + residue (tau = 1, 0.95), gain-aware closure.
import numpy as np, sys, time
from residue import closure_with_residue
from gac import gac
ids = [int(x) for x in sys.argv[1].split(",")]
for i in ids:
    t0 = time.time()
    Ws = list(np.load(f"W_off{i}.npy").astype(np.float64)); tr = np.load(f"truth_off{i}.npz"); mt = tr["m"][-1].astype(np.float64)
    r = closure_with_residue(Ws, tau=1.0); g = gac(Ws)
    mc = r[-1]["m_closure"]; gam1 = [d["gamma"] for d in r]; f = np.diff([0.0] + gam1)
    g95 = sum(f[l]*0.95**(len(f)-1-l) for l in range(len(f)))
    mse = lambda a: np.mean((a-mt)**2); c_or = mc @ (mc-mt)/(mc @ mc)
    sc = lambda v: 8*(mt @ (v-mt))/(mt @ mt)
    print(f"{i:4d} {mse(mc):.4e} {mse(mc*(1-gam1[-1]/8)):.4e} {mse(mc*(1-g95/8)):.4e} {mse(g[-1]['m']):.4e} {mse(mc*(1-c_or)):.4e}"
          f" | 8sc {sc(mc):+.4f} {sc(g[-1]['m']):+.4f} | gam {gam1[-1]:.4f} {g[-1]['gamma']:.4f} | noise {tr['avg_variance']/1e9:.1e} | {time.time()-t0:.0f}s", flush=True)
