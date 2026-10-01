"""Stage Q evaluator: run a numpy predictor on the shared benchmark sets and project to n = 1024.

Python use:
    import eval_q
    res = eval_q.evaluate(predict, ["w64_d16", "w128_d16", "w256_d16"], units_1024=120.0)
    # predict(W: np.ndarray (L, n, n) float32, x @ W convention) -> (L, n) per-neuron post-ReLU means
CLI:
    python eval_q.py --module my_design.py --func predict --sets w64_d16,w128_d16,w256_d16 --units 120 [--json out.json]
    python eval_q.py --baseline gauss --sets ...      # calibration baselines: gauss | mc
Reported per set: per-MLP final-layer MSE, raw = MSE - truth noise (avg_variance/N), all-layer MSE, per-layer MSE
(mean over MLPs), wall time per MLP.  Across sets: least-squares fit log(raw) = a - p log(n) over the d16 sets, the
extrapolated raw at n = 1024 (with the fit's 1-sigma band), and adjusted = raw_1024 * max(0.1, units/1024)
(1 unit = 2^31 FLOPs = one 1024^3 matmul, B = 2^41 = 1024 units).
"""
import argparse
import importlib.util
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bench  # noqa: E402


def eval_set(predict, name, mlps=None, verbose=True):
    S = bench.load_set(name)
    idx = range(len(S["seeds"])) if mlps is None else mlps
    rows = []
    for i in idx:
        W = bench.weights(S, i)
        t0 = time.time()
        p = np.asarray(predict(W), dtype=np.float64)
        dt = time.time() - t0
        T = S["means"][i]
        assert p.shape == T.shape, (p.shape, T.shape)
        lay = ((p - T) ** 2).mean(1)
        rows.append(dict(mlp=int(i), mse=float(lay[-1]), noise=float(S["noise"][i]), raw=float(lay[-1] - S["noise"][i]),
                         all_layers=float(lay.mean()), per_layer=lay.tolist(), finite=bool(np.isfinite(p).all()), wall=dt))
        if verbose:
            print(f"  {name} mlp {i}: mse {lay[-1]:.4e}  raw {lay[-1] - S['noise'][i]:.4e}  all-layer {lay.mean():.3e}  {dt:.1f}s",
                  flush=True)
    raw = np.array([r["raw"] for r in rows])
    out = dict(set=name, width=S["width"], depth=S["depth"], N=S["n_samples"], n_mlps=len(rows),
               mse=float(np.mean([r["mse"] for r in rows])), raw=float(raw.mean()),
               raw_se=float(raw.std(ddof=1) / np.sqrt(len(raw))) if len(raw) > 1 else float("nan"),
               noise=float(np.mean([r["noise"] for r in rows])),
               all_layers=float(np.mean([r["all_layers"] for r in rows])),
               per_layer=np.mean([r["per_layer"] for r in rows], 0).tolist(), per_mlp=rows)
    return out


def scaling_fit(results, target_n=1024, units_1024=None):
    pts = [(r["width"], r["raw"]) for r in results if r["depth"] == 16 and r["raw"] > 0 and r["width"] < target_n]
    fit = dict(points=pts)
    if len(pts) >= 2:
        x = np.log([p[0] for p in pts]); y = np.log([p[1] for p in pts])
        A = np.stack([np.ones_like(x), x], 1)
        coef, res, *_ = np.linalg.lstsq(A, y, rcond=None)
        a, b = coef
        pred = a + b * np.log(target_n)
        if len(pts) > 2:
            s2 = float(np.sum((y - A @ coef) ** 2) / (len(pts) - 2))
            cov = s2 * np.linalg.inv(A.T @ A)
            v = np.array([1.0, np.log(target_n)])
            sd = float(np.sqrt(v @ cov @ v))
        else:
            sd = float("nan")
        fit.update(p=float(-b), raw_1024=float(np.exp(pred)), raw_1024_band=[float(np.exp(pred - sd)), float(np.exp(pred + sd))])
        if units_1024 is not None:
            mult = max(0.1, units_1024 / 1024.0)
            fit.update(units_1024=units_1024, multiplier=mult, adjusted_1024=fit["raw_1024"] * mult)
    direct = [r for r in results if r["width"] == target_n and r["depth"] == 16]
    if direct:
        fit["raw_1024_measured"] = direct[0]["raw"]
        if units_1024 is not None:
            fit["adjusted_1024_measured"] = direct[0]["raw"] * max(0.1, units_1024 / 1024.0)
    return fit


def evaluate(predict, sets, units_1024=None, verbose=True):
    results = [eval_set(predict, s, verbose=verbose) for s in sets]
    fit = scaling_fit(results, units_1024=units_1024)
    if verbose:
        report(results, fit)
    return dict(sets=results, fit=fit)


def report(results, fit):
    print("| set | MLPs | N | final MSE | truth noise | raw | ± s.e. | all-layer MSE |")
    print("|---|---|---|---|---|---|---|---|")
    for r in results:
        print(f"| {r['set']} | {r['n_mlps']} | {r['N']:.0e} | {r['mse']:.3e} | {r['noise']:.1e} | {r['raw']:.3e} | "
              f"{r['raw_se']:.1e} | {r['all_layers']:.3e} |")
    if "p" in fit:
        print(f"width fit raw ∝ n^-{fit['p']:.2f}; extrapolated raw(1024) = {fit['raw_1024']:.3e} "
              f"(1σ band {fit['raw_1024_band'][0]:.2e}–{fit['raw_1024_band'][1]:.2e})")
    if "adjusted_1024" in fit:
        print(f"cost {fit['units_1024']} units = {fit['units_1024'] / 1024:.3f} B -> multiplier {fit['multiplier']:.3f}; "
              f"projected adjusted(1024) = {fit['adjusted_1024']:.3e}")
    if "raw_1024_measured" in fit:
        print(f"measured raw at 1024: {fit['raw_1024_measured']:.3e}")


# ---------------------------------------------------------------- calibration baselines (not designs)
def _phi(a):
    return np.exp(-0.5 * a * a) / np.sqrt(2 * np.pi)


def _Phi(a):
    from scipy.special import ndtr
    return ndtr(a)


def baseline_gauss(W):
    """Gaussian covariance closure (linearised cross-covariance, as in the starter kit): z ~ N(m, S) per layer,
    exact ReLU marginal mean/variance, post-ReLU covariance d(Phi) S_off d(Phi) + d(Var relu).  float64."""
    W = W.astype(np.float64)
    L, n, _ = W.shape
    out = []
    m = np.zeros(n); S = W[0].T @ W[0]
    for l in range(L):
        if l > 0:
            m = mu @ W[l]; S = W[l].T @ C @ W[l]
        v = np.maximum(np.diag(S), 1e-300); s = np.sqrt(v); a = m / s
        P = _Phi(a); p = _phi(a)
        mu = m * P + s * p
        sec = (m * m + v) * P + m * s * p
        out.append(mu)
        C = S * P[:, None] * P[None, :]
        np.fill_diagonal(C, np.maximum(sec - mu * mu, 0.0))
    return np.stack(out)


def baseline_mc(W, n_samples=6554, seed=0):
    """Plain Monte Carlo at the 0.1-floor budget of the n = 1024 competition (0.1 B / (2 L n^2) ≈ 6554 samples)."""
    rng = np.random.default_rng(seed)
    L, n, _ = W.shape
    h = rng.standard_normal((n_samples, n)).astype(np.float32)
    out = []
    for l in range(L):
        h = np.maximum(h @ W[l], 0.0)
        out.append(h.mean(0, dtype=np.float64))
    return np.stack(out)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--module"); ap.add_argument("--func", default="predict")
    ap.add_argument("--baseline", choices=["gauss", "mc"])
    ap.add_argument("--sets", required=True)
    ap.add_argument("--units", type=float)
    ap.add_argument("--json")
    a = ap.parse_args()
    if a.baseline:
        fn = dict(gauss=baseline_gauss, mc=baseline_mc)[a.baseline]
        units = dict(gauss=2 * 16 * 1.0, mc=102.4)[a.baseline] if a.units is None else a.units
    else:
        spec = importlib.util.spec_from_file_location("design", a.module)
        mod = importlib.util.module_from_spec(spec); sys.path.insert(0, os.path.dirname(os.path.abspath(a.module)))
        spec.loader.exec_module(mod)
        fn = getattr(mod, a.func); units = a.units
    res = evaluate(fn, a.sets.split(","), units_1024=units)
    if a.json:
        json.dump(res, open(a.json, "w"), indent=1)
