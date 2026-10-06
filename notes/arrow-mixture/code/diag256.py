import numpy as np, sys, time
sys.path.insert(0, "../num12")
from arrow import arrow
mse = lambda a, b: np.mean((a-b)**2)
Ws = list(np.load("../num12/W_n256_L16_s0.npy")); mt = np.load("../num12/truth_n256_L16_s0.npz")["m"]
ref = arrow(Ws, k=0)
print("closure %.3e" % mse(ref[-1], mt[-1]), flush=True)
for base, k, q, w, g in [("neuron",1,7,1,False), ("neuron",1,7,2,False), ("neuron",1,7,3,False), ("neuron",1,15,3,False), ("eig",1,7,3,False), ("rand",1,7,3,False), ("neuron",1,15,3,True)]:
    t0 = time.time(); a = arrow(Ws, k=k, q=q, w=w, base=base, gain=g)
    print(f"{base} k={k} q={q} w={w} gain={g}: {mse(a[-1], mt[-1]):.3e}  per-layer: " + " ".join(f"{mse(a[l], mt[l]):.1e}" for l in [3, 7, 11, 15]) + f"  ({time.time()-t0:.0f}s)", flush=True)
