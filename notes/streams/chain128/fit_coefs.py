"""Fit the per-layer coefficient table of the first-order closure (variant D) on atlases.

For each layer l the all-distinct kappa3(a_l) of the atlas, minus the Gaussian Hermite term to degree 4, is regressed
jointly over all given atlases on the seven diagram tensors of oracle_k3.residual_basis (built from the atlas's true
objects).  --space d21 fits instead in the transported space (the residual's contribution to D21(z_{l+1})).
    python fit_coefs.py --out coefs.npy [--space tensor|d21] A.npz B.npz ...
"""
import argparse
import numpy as np
import chain  # noqa: F401  (puts notes/experiments on the path)
import oracle_k3 as ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True); ap.add_argument("--space", default="tensor")
    ap.add_argument("atlases", nargs="+")
    a = ap.parse_args()
    zs = [np.load(p) for p in a.atlases]
    L = zs[0]["weights"].shape[0]
    coefs = np.zeros((L, 7)); coefs[:] = ok.CLOSURE_COEF
    for l in range(L - 1):
        Xs, ys = [], []
        for z in zs:
            o = ok.layer_objects(z, l)
            Kd = ok.all_distinct(o["K3a"])
            H = ok.hermite_model(o["C"], o["mu"], o["var"], 4)
            B = ok.residual_basis(o)
            if a.space == "d21":
                Wn = z["weights"][l + 1].astype(np.float64)
                ys.append(ok.offdiag(ok.transport_d21(Kd - H, Wn)).ravel())
                Xs.append(np.stack([ok.offdiag(ok.transport_d21(b, Wn)).ravel() for b in B], 1))
            else:
                ys.append((Kd - H).ravel()); Xs.append(np.stack([b.ravel() for b in B], 1))
        X = np.concatenate(Xs); y = np.concatenate(ys)
        live = np.abs(X).sum(0) > 0
        c = np.zeros(7)
        c[live], *_ = np.linalg.lstsq(X[:, live], y, rcond=None)
        if l == 0:
            c[~live] = 0.0
        coefs[l] = c
        r2 = 1 - np.sum((y - X @ c) ** 2) / np.sum(y ** 2)
        print(f"layer {l:2d}: R2 {r2:.3f}  coef " + " ".join(f"{v:+.2f}" for v in c), flush=True)
    np.save(a.out, coefs)


if __name__ == "__main__":
    main()
