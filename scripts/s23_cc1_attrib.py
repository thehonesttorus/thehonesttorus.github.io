"""Where CC1's remaining (incoherent) error lives (python scripts/s23_cc1_attrib.py KIK1 S21C DATA NET ...):
final RMS, coherent part, the part inherited through the penultimate means (Phi(alpha) W_L dmu_(L-1)), and the
readout's own residual; per-layer RMS errors; CC1 against the closure."""
import sys, os, glob, numpy as np
from scipy.special import ndtr
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from whest.cc1 import CC1
rms = lambda x: np.sqrt(np.mean(np.square(x)))
K1, SC, D = sys.argv[1], sys.argv[2], sys.argv[3]
for net in map(int, sys.argv[4:]):
    W = [np.asarray(w, dtype=np.float64) for w in np.load(f"{D}/W_off{net}.npy")]
    T = np.load(f"{D}/truth_off{net}.npz")["m"].astype(np.float64)
    k = np.load(f"{K1}/kikm_{net}.npz"); t, q = k["t"], k["q"]; al = t / np.sqrt(q - t * t)
    a = W[-1] @ (T[-2] / np.linalg.norm(T[-2])); x = a / a.std()
    line = [f"net {net}:"]
    for name, model in (("closure", CC1(collective=False)), ("cc1", CC1(Q=7))):
        out, _ = model.run(W); e = out[-1] - T[-1]; coh = np.polyval(np.polyfit(x, e, 8), x)
        inh = ndtr(al) * (W[-1] @ (out[-2] - T[-2]))
        line.append(f"  {name:8s} final {rms(e):.3e} (coherent {rms(coh):.2e}) | inherited {rms(inh):.2e}, readout residual {rms(e - inh):.2e} |"
                    f" layers 1/4/8/12/15: " + " ".join(f"{rms(out[l] - T[l]):.1e}" for l in (1, 4, 8, 12, 15)))
    print("\n".join(line), flush=True)
