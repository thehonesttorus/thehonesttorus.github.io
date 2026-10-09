"""Compare the chain's per-neuron third cumulants D3 with the Monte Carlo cache's kappa3(z), per layer, and measure
an oracle readout (MC kappa3 used at the readout only).   python scripts/diag_k3_mc.py DATA NET [opt=val ...]"""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from whest.k3chain import k3_chain
from whest.relu_gauss import phi, hermite_relu
D, net = sys.argv[1], int(sys.argv[2]); opts = {}
for a in sys.argv[3:]:
    k, v = a.split("="); opts[k] = v if k == "k4" else (float(v) if "." in v else int(v))
W = np.load(f"{D}/W_off{net}.npy"); mc = np.load(f"{D}/mccache_{net}.npz"); S = float(mc["S"])
mt = np.load(f"{D}/truth_off{net}.npz")["m"].astype(np.float64)
z1 = mc["s1"] / S; z2 = mc["s2"] / S; z3 = mc["s3"] / S; z4 = mc["s4"] / S
var_mc = z2 - z1 ** 2; k3_mc = z3 - 3 * z2 * z1 + 2 * z1 ** 3
k4_mc = z4 - 4 * z3 * z1 + 6 * z2 * z1 ** 2 - 3 * z1 ** 4 - 3 * var_mc ** 2
rec = {}; out, fl = k3_chain(W, opts, record=rec); L = W.shape[0]
print(f"opts {opts}: final MSE {np.mean((out[-1]-mt[-1])**2):.4e}  C/B {fl/2**41:.3f}")
print("layer | k3_mc rms | D3 rms | corr | slope | rel.err | MSE chain | MSE readout w/ MC k3 | w/ MC k3+k4 | w/ MC mean,var,k3,k4 (readout-only oracle)")
for l in range(L):
    r = rec[l]; k3 = k3_mc[l]; D3 = r["D3"]; var = r["var"]; mu = r["mu"]; sigma = np.sqrt(var); alpha = r["alpha"]
    ph = phi(alpha); d0 = hermite_relu(alpha, 0)[0]
    base = sigma * d0
    m_mc3 = base - k3 * alpha * ph / (6 * var)
    m_mc34 = m_mc3 + k4_mc[l] * (alpha ** 2 - 1) * ph / (24 * sigma ** 3)
    s2 = np.sqrt(var_mc[l]); a2 = z1[l] / s2; p2 = phi(a2); d02 = hermite_relu(a2, 0)[0]
    m_or = s2 * d02 - k3 * a2 * p2 / (6 * var_mc[l]) + k4_mc[l] * (a2 ** 2 - 1) * p2 / (24 * s2 ** 3)
    e = lambda v: np.mean((v - mt[l]) ** 2)
    print(f"{l:5d} | {np.sqrt(np.mean(k3**2)):8.2e} | {np.sqrt(np.mean(D3**2)):8.2e} | {np.corrcoef(D3,k3)[0,1]:6.3f} | {np.dot(D3,k3)/np.dot(k3,k3):6.3f} | {np.sqrt(np.mean((D3-k3)**2))/np.sqrt(np.mean(k3**2)):6.3f} | {e(out[l]):8.2e} | {e(m_mc3):8.2e} | {e(m_mc34):8.2e} | {e(m_or):8.2e}")
