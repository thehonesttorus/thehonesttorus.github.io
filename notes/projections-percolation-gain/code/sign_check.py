# Sign of the closure error (Prop. "the closure overestimates"): scale coefficient per layer and the fraction of
# significant per-neuron errors that are overestimates. Run with PYTHONPATH=../../trees-in-the-gaps/code
import numpy as np
from closure import closure
for n, L, s in [(256,4,0),(256,8,0),(256,16,0),(256,16,1),(256,16,2),(512,16,0),(1024,16,0)]:
    Ws = list(np.load(f"W_n{n}_L{L}_s{s}.npy")); tr = np.load(f"truth_n{n}_L{L}_s{s}.npz"); oc = closure(Ws)
    cs, fr = [], []
    for l in range(1, L):
        m = tr["m"][l]; mc = oc[l]["m"]; e = mc - m
        cs.append(mc @ e/(mc @ mc))
        sig = np.abs(e) > 3*np.sqrt(tr["v"][l]/tr["T"])
        fr.append(np.mean(e[sig] > 0) if sig.any() else np.nan)
    print(f"n={n} L={L} s={s}: c_l > 0 at {np.sum(np.array(cs)>0)}/{len(cs)} layers (min {min(cs):+.4f}); "
          f"final-layer overestimate fraction {fr[-1]:.2f}; min over layers {np.nanmin(fr):.2f}")
