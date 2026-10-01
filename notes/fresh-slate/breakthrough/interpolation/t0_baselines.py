"""T0: calibration at n=1024: Gaussian closure; closure x chi factor (radial mode); 'tree cover' = diagonal closure
(independent neurons, exact by additivity of cumulants up to O(n^-1/2) Edgeworth terms); annealed state evolution."""
import json, sys, time
import numpy as np
from common import bench, phi, Phi, chi_mean_ratio

def closure(W, diag_only=False, annealed_var=False):
    W = W.astype(np.float64); L, n, _ = W.shape
    out = []; m = np.zeros(n); S = W[0].T @ W[0]
    if diag_only: S = np.diag(np.diag(S))
    for l in range(L):
        if l > 0:
            m = mu @ W[l]
            if diag_only: S = np.diag((W[l] ** 2).T @ np.diag(C))
            else: S = W[l].T @ C @ W[l]
        v = np.maximum(np.diag(S), 1e-300)
        if annealed_var: v = np.full(n, (2.0 / n) * np.trace(C)) if l > 0 else v
        s = np.sqrt(v); a = m / s; P = Phi(a); p = phi(a)
        mu = m * P + s * p; sec = (m * m + v) * P + m * s * p
        out.append(mu)
        C = S * P[:, None] * P[None, :]
        np.fill_diagonal(C, np.maximum(sec - mu * mu, 0.0))
    return np.stack(out)

if __name__ == "__main__":
    name = sys.argv[1] if len(sys.argv) > 1 else "w1024_d16"
    S = bench.load_set(name); n = S["width"]; k = chi_mean_ratio(n)
    rows = []
    for i in range(len(S["seeds"])):
        W = bench.weights(S, i); T = S["means"][i]; nz = S["noise"][i]
        t0 = time.time()
        g = closure(W); d = closure(W, diag_only=True); an = closure(W, annealed_var=True)
        r = dict(mlp=i, noise=nz, gauss=((g[-1]-T[-1])**2).mean()-nz, gauss_chi=((k*g[-1]-T[-1])**2).mean()-nz,
                 tree=((d[-1]-T[-1])**2).mean()-nz, annealed_var=((an[-1]-T[-1])**2).mean()-nz,
                 truth_rms=float(np.sqrt((T[-1]**2).mean())), bias_gauss=float((g[-1]-T[-1]).mean()),
                 gauss_layers=((g-T)**2).mean(1).tolist(), wall=time.time()-t0)
        rows.append(r); print({kk: (f"{vv:.3e}" if isinstance(vv, float) else vv) for kk, vv in r.items() if kk != "gauss_layers"}, flush=True)
    for key in ["gauss", "gauss_chi", "tree", "annealed_var"]:
        print(key, f"{np.mean([r[key] for r in rows]):.3e}")
    json.dump(rows, open(f"results/t0_{name}.json", "w"), indent=1)
