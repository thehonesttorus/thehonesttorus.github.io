# Consistent against inconsistent all-oracle endpoint (V41_ORC_CONSIST), noise-free from the two Monte Carlo halves.
#   python endcmp.py NET DIR FREE.npy      (DIR holds tele_off<NET>_{A,B}_15.npy and teleinc_off<NET>_{A,B}_15.npy)
import sys, numpy as np
net, d, freep = int(sys.argv[1]), sys.argv[2], sys.argv[3]
mt = np.load(f"../official/truth_off{net}.npz")["m"][-1].astype(np.float64)
e = lambda p: np.load(p).astype(np.float64)[-1] - mt
ef = e(freep); E0 = float(ef @ ef)
for tag in ("tele", "teleinc"):
    a, b = e(f"{d}/{tag}_off{net}_A_15.npy"), e(f"{d}/{tag}_off{net}_B_15.npy")
    print(f"net {net} {tag:8s}: all-oracle output error, noise-free {100 * float(a @ b) / E0:6.2f}% of the free n MSE "
          f"(half A {100 * float(a @ a) / E0:.2f}%, half B {100 * float(b @ b) / E0:.2f}%)")
