"""'Not carried at all': can a chain drop the old pool and absorb its D21 effect into renormalised births?

For every layer l (target interface D21(z_l)), with sources in the published (no-AD) convention:
  share_w(l) = ||D21(Old_l^w)|| / ||D21(z_l)||,  Old^w = sources of age > w,  w = 1, 2, 4
and a least-squares replacement of D21(Old_l^w) by
  memoryless features : D21 transports W_l^{(x)3} of the birth diagrams of layer l-1 (oracle_k3: leading Gaussian Wick
                        rho^2 term, B1..B5 of residual_basis = D21(z)-hyperedge, D3, Gaussian rho^3, K22 diagrams;
                        B0 = the old content itself is excluded, B6 needs a --k4 atlas), plus their transposes
  + young features    : the D21 of each exactly carried young source (ages 1..w) and transposes
Coefficients are fitted per layer on one atlas and evaluated on others (held out); eps = ||resid|| / ||D21(z_l)||.

    python absorb.py FIT.npz [EVAL.npz ...]
"""
import sys, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "experiments"))
from tracker import load_atlas, run_sources, T3, d21  # noqa: E402
from oracle_k3 import residual_basis, hermite_model  # noqa: E402

WS = (1, 2, 4)


def collect(path):
    W, lay, N = load_atlas(path)
    L, n, _ = W.shape
    z = np.load(path)
    D, _, _ = run_sources(W, lay, ad=False)
    out = {}
    for l in range(2, L):
        den = np.linalg.norm(d21(np.asarray(lay[l]["K3z"], dtype=np.float64)))
        # birth diagrams of layer l-1 (needs the pre-activation objects of layer l-1 only)
        o = _objects(z, lay, l - 1)
        B = residual_basis(o)[1:]                      # drop B0 (= old content)
        B = [hermite_model(o["C"], o["mu"], o["var"], 4)] + B
        mem = []
        for b in B:
            f = d21(T3(b, W[l]))
            mem += [f, f.T]
        age = {l - s: D[s + 1, l] for s in range(0, l)}
        rec = dict(den=den, mem=mem, tgt={}, young={})
        for w in WS:
            rec["tgt"][w] = sum((age[a] for a in age if a > w), np.zeros((n, n)))
            yf = []
            for a in range(1, w + 1):
                if a in age:
                    yf += [age[a], age[a].T]
            rec["young"][w] = yf
        out[l] = rec
        print(f"  {path}: layer {l} collected", file=sys.stderr, flush=True)
    return out, n


def _objects(z, lay, l):
    """pre-activation objects of layer l for residual_basis (K22 from the pair moments; no K211)."""
    pre_s = z["pre_s"]
    mu, m2 = pre_s[0, l], pre_s[1, l]; var = m2 - mu ** 2
    M11 = z["pre_M11"][l].astype(np.float64); C = M11 - np.outer(mu, mu)
    M21 = z["pre_M21"][l].astype(np.float64); M22 = z["pre_M22"][l].astype(np.float64)
    Eu2u2 = (M22 - 2 * mu[None, :] * M21 - 2 * mu[:, None] * M21.T
             + np.outer(m2, mu ** 2) + np.outer(mu ** 2, m2) + 4 * np.outer(mu, mu) * M11 - 3 * np.outer(mu ** 2, mu ** 2))
    K22 = Eu2u2 - np.outer(var, var) - 2 * C ** 2; np.fill_diagonal(K22, 0.0)
    return dict(mu=mu, var=var, C=C, K3z=np.asarray(lay[l]["K3z"], dtype=np.float64), Phi=lay[l]["Phi"], K22=K22)


def fit(tgt, feats):
    X = np.stack([f.ravel() for f in feats], 1)
    c, *_ = np.linalg.lstsq(X, tgt.ravel(), rcond=None)
    return c


def ev(tgt, feats, c, den):
    X = np.stack([f.ravel() for f in feats], 1)
    return float(np.linalg.norm(tgt.ravel() - X @ c) / den)


def main():
    paths = sys.argv[1:]
    data = [collect(p) for p in paths]
    n = data[0][1]
    names = [os.path.basename(os.path.dirname(p)) for p in paths]
    print(f"width {n}; fit on {names[0]}, evaluated on {names}  (eps rel ||D21(z_l)||; 'share' = ||D21(Old^w)||/||D21(z_l)||)")
    for w in WS:
        print(f"\nw = {w}:  l | " + " | ".join(f"{nm}: share  mem  mem+young" for nm in names))
        F = data[0][0]
        for l in sorted(F):
            cm = fit(F[l]["tgt"][w], F[l]["mem"])
            cy = fit(F[l]["tgt"][w], F[l]["mem"] + F[l]["young"][w])
            row = []
            for (Dd, _) in data:
                r = Dd[l]
                row.append(f"{np.linalg.norm(r['tgt'][w])/r['den']:5.3f} {ev(r['tgt'][w], r['mem'], cm, r['den']):5.3f} "
                           f"{ev(r['tgt'][w], r['mem'] + r['young'][w], cy, r['den']):5.3f}")
            print(f"        {l:>2} | " + " | ".join(row), flush=True)


if __name__ == "__main__":
    main()
