"""T1: is the closure's coherent bias a fluctuating order parameter?  MC at n=1024 on bench networks.
Per layer l: tau_l(x) = |a_l(x)-mu_l|^2 / E|.|^2 (mu = bake truth). Records Var(tau), corr across layers, corr with |x|^2,
and, per neuron of layer l+1, the one-step readout defect on the empirical law:
  d_j = mean relu(z_j) - psi(mean z_j, var z_j)            (Gaussian readout with exact first two moments)
  g_j = mean relu(z_j) - E_tau psi(m_j + beta_j (tau-1), tau v_j')   (scale-mixture + coherent shift readout)"""
import sys, json, time
import numpy as np
from common import bench, phi, Phi

def psi(m, v):
    s = np.sqrt(np.maximum(v, 1e-300)); a = m / s
    return m * Phi(a) + s * phi(a)

name = sys.argv[1]; i = int(sys.argv[2]); N = int(sys.argv[3]); chunk = 4096
S = bench.load_set(name); W = bench.weights(S, i); T = S["means"][i]; L, n, _ = W.shape
rng = np.random.default_rng(123)
taus = []; zs = {}  # keep z of selected layers fully (float32) for the readout test
keep = [3, 7, 11, 15]
Z = {l: [] for l in keep}
t0 = time.time()
for c in range(N // chunk):
    h = rng.standard_normal((chunk, n)).astype(np.float32)
    x2 = (h.astype(np.float64) ** 2).sum(1)
    tl = [x2]
    for l in range(L):
        z = h @ W[l]
        if l in keep: Z[l].append(z)
        h = np.maximum(z, 0.0)
        tl.append(((h - T[l].astype(np.float32)) ** 2).sum(1, dtype=np.float64))
    taus.append(np.stack(tl, 1))
print("mc", time.time() - t0, flush=True)
tau = np.concatenate(taus); tau = tau / tau.mean(0)          # (N, L+1): col 0 = |x|^2, col l+1 = a_l
out = dict(var_tau=tau.var(0).tolist(), chi_ref=2.0 / n,
           corr_next=[float(np.corrcoef(tau[:, k], tau[:, k + 1])[0, 1]) for k in range(L)],
           corr_x=[float(np.corrcoef(tau[:, 0], tau[:, k])[0, 1]) for k in range(L + 1)])
print("Var tau (col0=|x|^2):", np.round(np.array(out["var_tau"]) * n / 2, 2), "(units of 2/n)")
print("corr(tau_l,tau_l+1):", np.round(out["corr_next"], 3))
print("corr(|x|^2,tau_l):", np.round(out["corr_x"], 3))
gh = np.polynomial.hermite_e.hermegauss(40)
for l in keep:
    z = np.concatenate(Z[l]).astype(np.float64); a = np.maximum(z, 0)
    tr = a.mean(0); m = z.mean(0); v = z.var(0)
    d = tr - psi(m, v)
    tp = tau[:, l]          # tau of a_{l-1} (col l = layer index l-1 output) -> mixing for z_l
    tc = tp - 1.0; vt = tc.var()
    beta = (tc @ (z - m)) / (len(tc) * vt)                    # coherent shift per neuron
    # residual variance after removing shift: v' with z = m + beta(tau-1) + sqrt(tau) xi, Var xi = v'
    vp = (v - beta ** 2 * vt) / 1.0
    # E over empirical tau via binning into quantiles (64 bins)
    qs = np.quantile(tp, (np.arange(64) + 0.5) / 64)
    g_pred = np.mean([psi(m + beta * (q - 1), q * vp) for q in qs], 0)
    gs_pred = np.mean([psi(m, q * v) for q in qs], 0)        # pure scale mixture, no shift
    dev_true = tr - T[l]   # MC noise of this run vs bake
    res = dict(layer=l, rms_defect_gauss=float(np.sqrt((d ** 2).mean())), mean_defect_gauss=float(d.mean()),
               rms_defect_gsm=float(np.sqrt(((tr - gs_pred) ** 2).mean())), mean_defect_gsm=float((tr - gs_pred).mean()),
               rms_defect_gsm_shift=float(np.sqrt(((tr - g_pred) ** 2).mean())), mean_defect_gsm_shift=float((tr - g_pred).mean()),
               rms_beta=float(np.sqrt((beta ** 2).mean())), rms_mc_vs_bake=float(np.sqrt((dev_true ** 2).mean())))
    print(res, flush=True); out[f"L{l}"] = res
json.dump(out, open(f"results/t1_{name}_{i}_N{N}.json", "w"), indent=1)
