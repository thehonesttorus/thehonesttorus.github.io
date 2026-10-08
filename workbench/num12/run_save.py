# run closure, TLP (J=1, rmax) and diag-TLP; save final-layer means
import numpy as np, sys, time
from closure import closure
from corrected3 import tree_closure
from diagsrc import diag_tlp
n, L, s, rmax = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
Ws = list(np.load(f"W_n{n}_L{L}_s{s}.npy"))
t0 = time.time(); oc = closure(Ws); t1 = time.time()
od = diag_tlp(Ws, rmax=2); t2 = time.time()
ot = tree_closure(Ws, rmax=rmax, J=1); t3 = time.time()
np.savez(f"est_n{n}_L{L}_s{s}.npz", closure=np.array([o["m"] for o in oc]), diag=np.array([o["m"] for o in od]), tlp=np.array([o["m"] for o in ot]))
print(n, L, s, "times closure %.1f diag %.1f tlp %.1f" % (t1-t0, t2-t1, t3-t2))
