"""Linear-response transfer coefficients at n = 1024: how a local error injected into the post-activation state of
layer l reaches the final-layer MSE, under the exact Gaussian closure dynamics (a proxy for the true linearised
dynamics; the transfer is a property of the propagation, not of the closure's own error).

Probes (central differences, +/- amplitude, so second-order terms cancel; dMSE = ||m_L(+) - m_L(-)||^2 / (4 n)):
  mean   : delta m_l iid Gaussian, rms 1                     -> K_mean[l]   = dMSE per (rms delta m)^2
  diag   : delta C(a_l)_ii iid Gaussian, rms 1                -> K_diag[l]   = dMSE per (rms delta var)^2
  off    : delta C(a_l) symmetric, zero diagonal, iid entries, Frobenius 1   -> K_off[l] = dMSE per ||dC||_F^2
  coh    : delta C(a_l) = C_off(a_l) / ||C_off(a_l)||_F (coherent shape)     -> K_coh[l] per ||dC||_F^2
  top    : delta C(a_l) = top-eigencomponent of C(a_l), normalised to Frobenius 1 -> K_top[l]
usage: python lr_probe.py MLP [probes]
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "bench"))
import bench  # noqa: E402
import gclose as g  # noqa: E402


def main(i, probes=("mean", "diag", "off", "coh", "top")):
    S = bench.load_set("w1024_d16")
    W = bench.weights(S, i).astype(np.float64)
    L, n, _ = W.shape
    base, states = g.run(W, "exact")
    rng = np.random.default_rng(77 + i)
    res = {p: [None] * L for p in probes}
    res["frob_Coff"] = [float(np.linalg.norm(C - np.diag(np.diag(C)))) for _, C in states]
    res["tr_C"] = [float(np.trace(C)) for _, C in states]
    res["top_eig"] = []
    for l in range(L - 1):
        m0, C0 = states[l]
        ev, evec = np.linalg.eigh(C0)
        res["top_eig"].append([float(x) for x in ev[-4:]])
        for p in probes:
            if p == "mean":
                dm = rng.standard_normal(n); dC = None; norm2 = 1.0
            elif p == "diag":
                dm = None; dC = np.diag(rng.standard_normal(n)); norm2 = 1.0
            elif p == "off":
                X = rng.standard_normal((n, n)); X = (X + X.T) / np.sqrt(2); np.fill_diagonal(X, 0)
                dm = None; dC = X / np.linalg.norm(X); norm2 = 1.0
            elif p == "coh":
                X = C0 - np.diag(np.diag(C0)); dm = None; dC = X / np.linalg.norm(X); norm2 = 1.0
            elif p == "top":
                v = evec[:, -1]; X = np.outer(v, v); dm = None; dC = X / np.linalg.norm(X); norm2 = 1.0
            # amplitude: small relative to the state
            amp = 1e-3 if p in ("mean", "diag") else 1e-2
            outs = []
            for sgn in (+1, -1):
                m1 = m0 + (sgn * amp * dm if dm is not None else 0)
                C1 = C0 + (sgn * amp * dC if dC is not None else 0)
                if l == L - 1:
                    outs.append(m1)
                else:
                    o, _ = g.run(W, "exact", start=(l + 1, m1, C1))
                    outs.append(o[-1])
            d = (outs[0] - outs[1]) / (2 * amp)
            if p == "mean":
                norm2 = float(np.mean(dm ** 2))
            elif p == "diag":
                norm2 = float(np.mean(np.diag(dC) ** 2))
            res[p][l] = float(np.mean(d ** 2) / norm2)
            print(f"mlp {i} l {l} {p}: K = {res[p][l]:.3e}", flush=True)
    json.dump(res, open(os.path.join(HERE, "results", f"lr_mlp{i}.json"), "w"), indent=1)


if __name__ == "__main__":
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    i = int(sys.argv[1])
    pr = sys.argv[2].split(",") if len(sys.argv) > 2 else ("mean", "diag", "off", "coh", "top")
    main(i, pr)
