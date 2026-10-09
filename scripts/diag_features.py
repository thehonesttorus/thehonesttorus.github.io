"""Is the per-neuron error of a chain predictable from the chain's own cheap per-neuron features? Fit a ridge regression
of the final-layer error on O(n)-per-layer features from the last two layers on one network, test on the other.
  python scripts/diag_features.py DATA [chain opts...]   (opts as in run_kprop3c.py; window=99 is the exact chain)"""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from whest.kprop3c import kprop3c_chain
from whest.kprop3 import wick
D = sys.argv[1]; opts = {}
for a in sys.argv[2:]:
    k, v = a.split("="); opts[k] = v if k in ("oldmode", "tier", "specres") else (float(v) if "." in v else int(v))
def feats(rec, l):
    r = rec[l]; m, var = r["m"], r["var"]; sd = np.sqrt(var); al = m / sd
    F = [m, sd, al, wick(m, var, 1, 1), wick(m, var, 2, 1), wick(m, var, 3, 1), r["K3_3"], r["K3h_3"], r["mu_h"], np.diag(r["Ch"]),
         np.sqrt(np.sum(r["K3_21"] ** 2, 1)), np.sqrt(np.sum(r["K3_21"] ** 2, 0)), np.sqrt(np.sum(r["K3h_21"] ** 2, 1)),
         np.full(len(m), r["c4"]), r["K3_3"] / sd ** 3, r["K3_21"].sum(1), r["K3_21"].sum(0)]
    return np.stack(F, 1)
data = {}
for net in (0, 1):
    W = np.load(f"{D}/W_off{net}.npy").astype(np.float64); L = W.shape[0]
    truth = np.load(f"{D}/truth_off{net}.npz")["m"].astype(np.float64)
    rec = {}; out = kprop3c_chain(W, dict(opts), record=rec); err = out[-1] - truth[-1]
    X = np.concatenate([feats(rec, L - 1), feats(rec, L - 2)], 1); data[net] = (X, err)
    print(f"net {net}: raw MSE {np.mean(err**2):.4e}; collective share of the error {np.dot(err, rec[L-1]['mu_h'])**2/np.dot(err,err)/np.dot(rec[L-1]['mu_h'],rec[L-1]['mu_h']):.3f}")
def design(X, deg):
    Xs = (X - mu) / sdv; cols = [np.ones(len(X))] + [Xs[:, j] for j in range(Xs.shape[1])]
    if deg >= 2: cols += [Xs[:, i] * Xs[:, j] for i in range(Xs.shape[1]) for j in range(i, Xs.shape[1])]
    return np.stack(cols, 1)
for deg in (1, 2):
    for tr, te in ((0, 1), (1, 0)):
        Xtr, etr = data[tr]; Xte, ete = data[te]; mu, sdv = Xtr.mean(0), Xtr.std(0) + 1e-12
        A, B = design(Xtr, deg), design(Xte, deg)
        for lam in (1e-3, 1e-1, 10.0):
            beta = np.linalg.solve(A.T @ A + lam * np.eye(A.shape[1]), A.T @ etr)
            print(f"deg {deg} ridge {lam:g}: fit on net {tr} -> net {te}: MSE {np.mean((ete - B @ beta)**2):.4e} (train {np.mean((etr - A @ beta)**2):.4e}); raw {np.mean(ete**2):.4e}")
