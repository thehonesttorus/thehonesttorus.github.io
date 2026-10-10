"""Tier-2 oracle: the Gaussian closure with its collective sector corrected at every layer by three Monte Carlo scalars
(python scripts/s21_tier2_oracle.py RESULTS DATA NET ...).

After each layer's closure step (m, K) of the post-activation, with mu_hat the mean direction:
  a: shift m along mu_hat so that E c = mu_hat . m matches;
  b: + rank-one update of K along mu_hat so that Var c = mu_hat K mu_hat matches (the collective variance);
  c: + rescale the perp block P K P so that E r^2 = Tr(P S^u) matches.
Everything else (per-neuron, incoherent) is the closure's own. Reports the final error and its coherent part, and the
closure's own collective variance error per layer."""
import sys, os, glob, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from whest.relu_gauss import relu_moments
from scripts.s21_collective import MIX
rms = lambda x: np.sqrt(np.mean(np.square(x)))
R, D = sys.argv[1], sys.argv[2]
for net in map(int, sys.argv[3:]):
    fs = sorted(glob.glob(f"{R}/s21c_{net}_*.npz")); N = 0; mix = 0
    for f in fs: z = np.load(f); N += int(z["N"]); mix = mix + z["mix"]
    M = {ij: mix[:, k] / N for k, ij in enumerate(MIX)}
    Ec, Vc, Er2 = M[(1, 0)], M[(2, 0)] - M[(1, 0)] ** 2, M[(0, 1)]
    truth = np.load(f"{D}/truth_off{net}.npz")["m"].astype(np.float64); L, n = truth.shape
    W = [np.asarray(w, dtype=np.float64) for w in np.load(f"{D}/W_off{net}.npy")]
    a = W[-1] @ (truth[-2] / np.linalg.norm(truth[-2])); x = a / a.std()
    res = {}
    for mode in ("none", "a", "ab", "abc"):
        mu_h = np.zeros(n); C_h = np.eye(n); verr = []
        for l in range(L):
            mu = W[l] @ mu_h; C = W[l] @ C_h @ W[l].T; m, K, _, _ = relu_moments(mu, C, 10)
            u = m / np.linalg.norm(m); verr.append(u @ K @ u / Vc[l] - 1)
            if l < L - 1 and mode != "none":
                m = m + (Ec[l] - u @ m) * u
                if "b" in mode:
                    Ku = K @ u; K = K + (Vc[l] - u @ Ku) * np.outer(u, u)
                if "c" in mode:
                    Pm = m - (u @ m) * u; KP = K - np.outer(K @ u, u) - np.outer(u, u @ K) + (u @ K @ u) * np.outer(u, u)
                    lam = (Er2[l] - Pm @ Pm) / np.trace(KP); K = K + (lam - 1) * KP
            mu_h, C_h = m, K
        e = mu_h - truth[-1]; res[mode] = (rms(e), rms(np.polyval(np.polyfit(x, e, 8), x)), verr)
    print(f"net {net}: final err (coherent part): " + " | ".join(f"{k}: {v[0]:.2e} ({v[1]:.2e})" for k, v in res.items()))
    print("   closure's own collective-variance error mu_hat K mu_hat / Var c - 1 by layer: " + " ".join(f"{v:+.3f}" for v in res["none"][2]), flush=True)
