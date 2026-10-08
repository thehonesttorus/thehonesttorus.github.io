# Local (one-step) closure error per layer from the TRUE pre-activation moments, its scale component, and how much the
# scale-mixture shape (g/8) sigma phi(a) (1+a^2) with the measured kappa3 coefficient g_l explains; Edgeworth kappa3 and
# "kappa4 = 3 g sigma^4" terms.  Also the closure's own (from-scratch) output scale error and GAC's.
import numpy as np, sys
from scipy.special import ndtr
sys.path.insert(0, "../num12")
from closure import relu_coeffs, closure
from gac import gac
from pairvar import y_cumulants, phi
net = int(sys.argv[1]) if len(sys.argv) > 1 else 0
D = np.load(f"mc_cum_off{net}.npz"); L, n = D["s1y"].shape
Wc = np.load(f"../official/W_off{net}.npy").astype(np.float64); mt = np.load(f"../official/truth_off{net}.npz")["m"].astype(np.float64)
gm = np.load(f"k3channel_off{net}.npz")["gm"] if len(sys.argv) < 3 else None
print(f"net {net}: local closure error per layer (closure step from the TRUE mu, sigma of z_l minus true mean of h_l)")
print(" l | local MSE | scale c_l (8c) | MSE explained: kappa3 Edgeworth / + kappa4=3 g s^4 / scale-mixture shape (g/8) s phi (1+a^2) | 8*scale of each")
cs = []
for l in range(L):
    Y = y_cumulants(D, l); mu, S, k3 = Y["mu"], Y["S"], Y["k3"]; sig = np.sqrt(np.diag(S)); a = mu / sig
    mG = sig * (a * ndtr(a) + phi(a)); err = mG - mt[l]; m = mt[l]
    sc = lambda v: (m @ v) / (m @ m); ex = lambda v: 1 - np.mean((err - v)**2) / np.mean(err**2)
    A = relu_coeffs(mu, sig, 6)
    e3 = -(k3 / 6) * A[3] / sig**3
    g = gm[l] if gm is not None else 0.0
    e4 = -(3 * g * sig**4 / 24) * A[4] / sig**4
    esm = (g / 8) * sig * phi(a) * (1 + a * a)
    cs.append(sc(err))
    print(f" {l:2d} | {np.mean(err**2):.2e} | {sc(err):+.6f} ({8*sc(err):+.5f}) | {ex(e3):.3f} / {ex(e3+e4):.3f} / {ex(esm):.3f} | {8*sc(e3):+.5f} {8*sc(e3+e4):+.5f} {8*sc(esm):+.5f}")
oc = closure(list(Wc)); gg = gac(list(Wc))
mL = mt[-1]; scL = lambda v: (mL @ (v - mL)) / (mL @ mL)
print(f"\noutput: closure from scratch: MSE {np.mean((oc[-1]['m']-mL)**2):.3e}, 8*scale {8*scL(oc[-1]['m']):+.5f};  GAC: MSE {np.mean((gg[-1]['m']-mL)**2):.3e}, 8*scale {8*scL(gg[-1]['m']):+.5f}, 8(1-E[G]) {8*(1-gg[-1]['g1']):.5f}")
print(f"sum of local scale errors (critical transport): 8*sum c_l = {8*np.sum(cs):+.5f}; with 0.9/layer retention: {8*np.sum([c*0.9**(L-1-l) for l, c in enumerate(cs)]):+.5f}")
