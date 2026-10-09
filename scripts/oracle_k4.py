"""Oracle: inject Monte Carlo per-neuron kappa4(z_l) (from the cache) into chain v2's readouts at every layer.
  python scripts/oracle_k4.py DATA NET [k4var=0/1]"""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from whest.k3chain2 import k3_chain2
D, net = sys.argv[1], int(sys.argv[2]); opts = {}
for a in sys.argv[3:]:
    k, v = a.split("="); opts[k] = int(v)
mc = np.load(f"{D}/mccache_{net}.npz"); S = float(mc["S"]); mt = np.load(f"{D}/truth_off{net}.npz")["m"].astype(np.float64)
z1 = mc["s1"] / S; z2 = mc["s2"] / S; z3 = mc["s3"] / S; z4 = mc["s4"] / S; var = z2 - z1 ** 2
k4 = z4 - 4 * z3 * z1 + 6 * z2 * z1 ** 2 - 3 * z1 ** 4 - 3 * var ** 2
k3 = z3 - 3 * z2 * z1 + 2 * z1 ** 3
W = np.load(f"{D}/W_off{net}.npy")
for name, o in [("base", dict(opts)), ("k4 oracle (readouts, all layers)", dict(opts, k4diag=k4)), ("k4 oracle last layer only", dict(opts, k4diag=np.vstack([np.zeros_like(k4[:-1]), k4[-1:]])))]:
    out, fl = k3_chain2(W, o); per = np.mean((out - mt) ** 2, axis=1)
    print(f"{name:40s}: final {per[-1]:.3e} | " + " ".join(f"{p:.1e}" for p in per), flush=True)
print("MC kappa4 rms per layer:", " ".join(f"{x:.2e}" for x in np.sqrt(np.mean(k4 ** 2, 1))), " noise ~", f"{np.sqrt(96 * np.mean(var[-1] ** 4) / S):.1e}")
print("excess kurtosis rms per layer:", " ".join(f"{x:.3f}" for x in np.sqrt(np.mean((k4 / var ** 2) ** 2, 1))))
