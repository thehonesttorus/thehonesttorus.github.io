# Localization diagnostics from the chain dump (chaindump_{n}.npz) and the MC cache (mccache_{n}.npz):
#   python ana_loc.py DIR NET [NET ...]
# E1 error vs gate margin alpha = mu / sigma (per layer: share of the squared mean error carried by |alpha| bins)
# E2 the chain's pre-activation mean / variance error against Monte Carlo (statistic-level accuracy by layer)
# E3 first-chaos share of each neuron's post-activation variance, |G[:, i]|^2 / Var(x_i), by layer
# E4 spectrum of the pre-activation covariance (MC): participation ratio and the share of the top 16 / 128 modes
# E5 injected error per layer: e_l - diag(Phi(alpha_l)) W_l e_(l-1) (first-chaos transport of the previous error)
#   python ana_loc.py --save DIR NET   writes $OUT/analoc_{NET}.npz;   python ana_loc.py --agg FILE...   aggregates
import os, sys, numpy as np
from math import erf, sqrt
SAVE = sys.argv[1] == "--save"; AGG = sys.argv[1] == "--agg"
if SAVE or AGG:
    sys.argv.pop(1)
D = sys.argv[1]; nets = [] if AGG else [int(x) for x in sys.argv[2:]]
Phi = np.vectorize(lambda z: 0.5 * (1 + erf(z / sqrt(2))))
bins = [0, 0.5, 1, 1.5, 2, 3, 99]
agg = {k: [] for k in ("e1", "e2m", "e2v", "e3", "e4", "e5", "mse")}
for net in nets:
    cd = np.load(f"{D}/chaindump_{net}.npz"); mc = np.load(f"{D}/mccache_{net}.npz")
    mt = np.load(f"{D}/truth_off{net}.npz")["m"].astype(np.float64)
    W = np.load(f"{D}/W_off{net}.npy").astype(np.float64)       # (L, n_out, n_in): h_l = W[l] @ x_l
    S = float(mc["S"]); L, n = mt.shape
    out = cd["out"]; err = out - mt
    m1 = mc["s1"] / S; m2 = mc["s2"] / S
    var_mc = m2 - m1 ** 2
    e1 = np.zeros((L, len(bins) - 1)); e2m = np.zeros(L); e2v = np.zeros(L); e3 = np.zeros(L); e4 = np.zeros((L, 3)); e5 = np.zeros(L)
    xm = mc["xs1"] / S
    for l in range(L):
        al = np.abs(cd["mu"][l]) / np.sqrt(np.maximum(cd["var"][l], 1e-12))
        tot = np.sum(err[l] ** 2)
        for b in range(len(bins) - 1):
            sel = (al >= bins[b]) & (al < bins[b + 1])
            e1[l, b] = np.sum(err[l][sel] ** 2) / tot
        e2m[l] = np.sqrt(np.mean((cd["mu"][l] - m1[l]) ** 2)); e2v[l] = np.sqrt(np.mean((cd["var"][l] / var_mc[l] - 1) ** 2))
        G = mc["G"][l].astype(np.float64)                          # E[X_j x_(l+1),i]  (j input, i neuron)
        # Var of the post-activation from the uncentred second moment diag(E[x x^T]) is not cached; use E[h^2] split:
        # first-chaos share of the PRE-activation h_l: |E[grad h]|^2 / Var(h); E[grad h_l] = W_l G_(l-1)
        Gin = np.eye(n) if l == 0 else mc["G"][l - 1].astype(np.float64)
        gh = Gin @ W[l].T                                           # (input j, neuron i) = E[X_j h_l,i]
        e3[l] = np.median(np.sum(gh ** 2, 0) / np.maximum(var_mc[l], 1e-12))
        Hc = mc["Hc"][l].astype(np.float64) - np.outer(m1[l], m1[l])
        ev = np.clip(np.linalg.eigvalsh(Hc)[::-1], 0, None)
        e4[l] = [ev.sum() ** 2 / np.sum(ev ** 2) / n, ev[:16].sum() / ev.sum(), ev[:128].sum() / ev.sum()]
        if l > 0:
            pg = Phi(cd["mu"][l] / np.sqrt(np.maximum(cd["var"][l], 1e-12)))
            inj = err[l] - pg * (W[l] @ err[l - 1])
            e5[l] = np.mean(inj ** 2) / np.mean(err[l] ** 2)
    for k, v in (("e1", e1), ("e2m", e2m), ("e2v", e2v), ("e3", e3), ("e4", e4), ("e5", e5), ("mse", np.mean(err ** 2, 1))):
        agg[k].append(v)
    if SAVE:
        np.savez(f"{os.environ.get('OUT', '.')}/analoc_{net}.npz", **{k: v[-1] for k, v in agg.items()})
    print(f"net {net}: final MSE {np.mean(err[-1] ** 2):.3e}; MC-vs-truth mean check (layer 15) {np.mean((xm[-1] - mt[-1]) ** 2):.1e}", flush=True)
if AGG:
    for f in sys.argv[1:]:
        z = np.load(f)
        for k in agg:
            agg[k].append(z[k])
    print(f"{len(agg['mse'])} networks")
A = {k: np.mean(v, 0) for k, v in agg.items()}
print("\nlayer  MSE      | share of sq. error by |alpha| bins " + " ".join(f"[{bins[b]},{bins[b+1]})" for b in range(len(bins) - 1))
      + " | rms(mu err)  rms(var rel err) | 1st-chaos share | PR/n top16 top128 | injected/total")
for l in range(len(A["mse"])):
    print(f"{l:2d}  {A['mse'][l]:.2e} | " + " ".join(f"{x:6.3f}" for x in A["e1"][l])
          + f" | {A['e2m'][l]:.2e}  {A['e2v'][l]:.2e} | {A['e3'][l]:.3f} | " + " ".join(f"{x:.3f}" for x in A["e4"][l])
          + f" | {A['e5'][l]:.3f}")
