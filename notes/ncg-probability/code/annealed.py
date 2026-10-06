# (E7) Is the chain's one-step residual annealed (a function of the per-neuron state, recoverable at O(n^2)) or quenched?
# Residual r_l(i) = m*_l(i) - chain_l(i) with the true mean injected at every layer (dump_oracleall).  Features per neuron:
# functions of the chain's own per-neuron state (mu, sigma, a, phi, Phi, D3, |D21| row norm, correlation row energy,
# previous true mean), polynomial degree 2 in the standardised features.  Ridge fit on one network, tested on the other
# and across layers; the out-of-sample R^2 is the annealed share.
import pickle, numpy as np, itertools
from scipy.stats import norm
def load(net):
    D = {d["layer"]: d for d in pickle.load(open(f"../k3work/dump_oracleall_off{net}.pkl", "rb"))}
    pred = np.load(f"../k3work/pred_oracleall_off{net}.npy"); mt = np.load(f"../official/truth_off{net}.npz")["m"].astype(np.float64)
    X, Y, lab = [], [], []
    for l in sorted(D):
        d = D[l]
        if d["C_pre"] is None or l < 3: continue
        mu, var = d["mu"], d["var"]; s = np.sqrt(var); a = mu / s; ph, Ph = norm.pdf(a), norm.cdf(a)
        C = d["C_pre"]; R = C / np.outer(s, s); np.fill_diagonal(R, 0.0); srow = (R**2).sum(1)
        D3 = d["D3"]; D21n = np.abs(d["D21"]).sum(1)
        feats = np.column_stack([mu, s, a, ph, Ph, D3 / s**3, D21n / s**3, srow, mt[l - 1], ph * s, a * ph, (a**2 - 1) * ph, Ph * s])
        X.append(feats); Y.append(mt[l] - pred[l]); lab.append(np.full(len(mu), l))
    return np.vstack(X), np.concatenate(Y), np.concatenate(lab)
def design(X, mean, std, deg2=True):
    Z = (X - mean) / std
    cols = [np.ones(len(Z))] + [Z[:, j] for j in range(Z.shape[1])]
    if deg2:
        cols += [Z[:, i] * Z[:, j] for i, j in itertools.combinations_with_replacement(range(Z.shape[1]), 2)]
    return np.column_stack(cols)
def ridge(A, y, lam):
    return np.linalg.solve(A.T @ A + lam * np.eye(A.shape[1]), A.T @ y)
r2 = lambda y, yh: 1 - np.mean((y - yh)**2) / np.mean(y**2)
X0, Y0, L0 = load(0); X1, Y1, L1 = load(1)
mean, std = X0.mean(0), X0.std(0) + 1e-12
for deg2 in (False, True):
    A0, A1 = design(X0, mean, std, deg2), design(X1, mean, std, deg2)
    for lam in (1e-2, 1e0, 1e2):
        b = ridge(A0, Y0, lam)
        print(f"deg2={deg2} lam={lam:g}: fit net0 -> in-sample R2 {r2(Y0, A0 @ b):.4f}, out-of-sample net1 R2 {r2(Y1, A1 @ b):.4f}; "
              f"reverse: {r2(Y1, A1 @ ridge(A1, Y1, lam)):.4f} / {r2(Y0, A0 @ ridge(A1, Y1, lam)):.4f}", flush=True)
# per-layer, same network, leave-one-layer-out (fit on other layers of net 0, test on layer l of net 0)
A0 = design(X0, mean, std, True)
print("leave-one-layer-out on net 0 (deg2, lam=1):")
for l in sorted(set(L0)):
    tr = L0 != l; b = ridge(A0[tr], Y0[tr], 1.0)
    print(f"  layer {l:2d}: out-of-layer R2 {r2(Y0[~tr], A0[~tr] @ b):+.4f}   (rms r {np.sqrt(np.mean(Y0[~tr]**2)):.2e})")
