import numpy as np, time
from closure import closure
from residue import closure_with_residue
from corrected3 import tree_closure
Ws = list(np.load("W_n1024_L16_s0.npy"))
t0 = time.time(); oc = closure(Ws); t1 = time.time()
r95 = closure_with_residue(Ws, tau=0.95); r100 = closure_with_residue(Ws, tau=1.0); t2 = time.time()
print("closure %.1fs residue %.1fs" % (t1-t0, (t2-t1)/2), flush=True)
np.savez("est_n1024_L16_s0_fast.npz", closure=np.array([o["m"] for o in oc]), res95=np.array([o["m"] for o in r95]),
         res100=np.array([o["m"] for o in r100]), gamma=np.array([o["gamma"] for o in r95]))
t3 = time.time(); ot = tree_closure(Ws, rmax=4, J=1); t4 = time.time()
print("tlp %.1fs" % (t4-t3), flush=True)
np.savez("est_n1024_L16_s0_tlp.npz", tlp=np.array([o["m"] for o in ot]))
