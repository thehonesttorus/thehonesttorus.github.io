# Official networks: GAC + fresh kappa_3 trees (D3, P3 j<=2, T3 J=1).
import numpy as np, sys, time
from gac_k3 import gac_k3
for i in [int(x) for x in sys.argv[1].split(",")]:
    t0 = time.time(); Ws = list(np.load(f"W_off{i}.npy").astype(np.float64)); mt = np.load(f"truth_off{i}.npz")["m"][-1].astype(np.float64)
    g = gac_k3(Ws, jmax=2)[-1]["m"]
    print(f"{i:4d} {np.mean((g-mt)**2):.4e} {8*(mt @ (g-mt))/(mt @ mt):+.4f} {time.time()-t0:.0f}s", flush=True)
