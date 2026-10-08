# Ceiling of the tree family on an official network: TLP (note IX) with linear-response transport, plus oracle scale.
import numpy as np, sys, time
from corrected3 import tree_closure
tag = sys.argv[1]; rmax = int(sys.argv[2])
Ws = list(np.load(f"W_{tag}.npy").astype(np.float64)); mt = np.load(f"truth_{tag}.npz")["m"][-1].astype(np.float64)
t0 = time.time(); m = tree_closure(Ws, rmax=rmax, J=3)[-1]["m"]
c = m @ (m - mt)/(m @ m)
print(f"{tag} TLP rmax={rmax}: MSE {np.mean((m-mt)**2):.3e}, with oracle scale {np.mean((m*(1-c)-mt)**2):.3e} (8c={8*c:+.4f}) [{time.time()-t0:.0f}s]", flush=True)
