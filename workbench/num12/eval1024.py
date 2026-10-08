import numpy as np
tr = np.load("truth_n1024_L16_s0.npz"); m = tr["m"][-1]; v = tr["v"][-1]; T = int(tr["T"])
f = np.load("est_n1024_L16_s0_fast.npz"); t = np.load("est_n1024_L16_s0_tlp.npz")
mc = f["closure"][-1]; mse = lambda x: np.mean((x-m)**2)
c = mc @ (mc-m)/(mc @ mc)
tl = t["tlp"][-1]; ct = tl @ (tl-m)/(tl @ tl)
print(f"MC@B {v.mean()/65536:.3e}  noise {v.mean()/T:.2e}")
print(f"closure {mse(mc):.3e}  res tau=1 {mse(f['res100'][-1]):.3e}  res tau=.95 {mse(f['res95'][-1]):.3e}  oracle scale {mse(mc*(1-c)):.3e} (c={c:.4f}, gamma/8 tau.95 = {f['gamma'][-1]/8:.4f})")
print(f"TLP {mse(tl):.3e}  TLP+oracle scale {mse(tl*(1-ct)):.3e} (c={ct:.4f})")
for l in [1, 3, 7, 11, 15]:
    ml = tr["m"][l]; print(f"  layer {l+1}: closure {np.mean((f['closure'][l]-ml)**2):.2e}  res.95 {np.mean((f['res95'][l]-ml)**2):.2e}  TLP {np.mean((t['tlp'][l]-ml)**2):.2e}  MC@B {tr['v'][l].mean()/65536:.2e}")
