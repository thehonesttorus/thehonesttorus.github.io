"""Attribute the Gaussian closure's final-mean error to the two tiers (python scripts/s21_tier_attrib.py RESULTS DATA NET ...).

Coherent tier: Ybar(a) = E_x G(a c, s^2 r^2) depends only on the law of (c, r^2) at the penultimate layer.
  (i)  the closure's collective scalars E c = mu_hat . m_hat and E r^2 = Tr S^u_hat - E c^2... against Monte Carlo;
  (ii) coherent-tier replacement: closure output + [Ybar_true(a) - fit of the closure output on a] removes the coherent
       error exactly; the remainder is the closure's per-neuron (incoherent) error.
  (iii) two-scalar repair: shift the closure's penultimate mean along mu_hat and rescale its covariance so that E c and
       E||h||^2 match Monte Carlo, then re-read the final layer: how much of the coherent error do two numbers carry."""
import sys, os, glob, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from whest import gate_channel as gc
from whest.relu_gauss import relu_moments
from scripts.s21_collective import MIX
rms = lambda x: np.sqrt(np.mean(np.square(x)))
R, D = sys.argv[1], sys.argv[2]
for net in map(int, sys.argv[3:]):
    fs = sorted(glob.glob(f"{R}/s21c_{net}_*.npz")); N = 0; mix = 0; yb = 0
    for f in fs:
        z = np.load(f); N += int(z["N"]); mix = mix + z["mix"]; yb = yb + z["yb_15"]
    M = {ij: mix[:, k] / N for k, ij in enumerate(MIX)}; Ytrue = yb / N; a = np.load(fs[0])["a_15"]; x = a / a.std()
    truth = np.load(f"{D}/truth_off{net}.npz")["m"].astype(np.float64); L, n = truth.shape
    W = [np.asarray(w, dtype=np.float64) for w in np.load(f"{D}/W_off{net}.npy")]
    st = gc.closure_states(W); pen = st[L - 2]; mh = truth[L - 2] / np.linalg.norm(truth[L - 2])
    m_hat, Kh = pen["m"], pen["Kh"]; S_hat = Kh + np.outer(m_hat, m_hat)
    Ec_cl, Eh2_cl = mh @ m_hat, np.trace(S_hat); Ec_mc, Eh2_mc = M[(1, 0)][L - 2], M[(2, 0)][L - 2] + M[(0, 1)][L - 2]
    y_cl = st[-1]["m"]; e = y_cl - truth[-1]; fit = np.polyval(np.polyfit(x, y_cl, 8), x)
    y_rep = y_cl + (Ytrue - fit)
    # (iii) two-scalar repair at the penultimate layer
    m2 = m_hat + (Ec_mc - Ec_cl) * mh
    lam = (Eh2_mc - m2 @ m2) / np.trace(Kh); Kh2 = Kh * lam
    mu_z = W[-1] @ m2; C_z = W[-1] @ Kh2 @ W[-1].T; y2, _, _, _ = relu_moments(mu_z, C_z, 10)
    e2 = y2 - truth[-1]; coh = lambda r: rms(np.polyval(np.polyfit(x, r, 8), x))
    print(f"net {net}: closure E c rel err {Ec_cl / Ec_mc - 1:+.2e}, E|h|^2 rel err {Eh2_cl / Eh2_mc - 1:+.2e} | final err {rms(e):.2e} (coherent {coh(e):.2e}) |"
          f" coherent tier replaced: {rms(y_rep - truth[-1]):.2e} | two-scalar repair: {rms(e2):.2e} (coherent {coh(e2):.2e})", flush=True)
