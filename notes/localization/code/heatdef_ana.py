# Analysis of the heat-defect experiment (note XLII):   python heatdef_ana.py RESULTDIR TRUTHDIR H NET [NET ...]
# delta = 2 tr D, estimated as 2 (n/K) sum_r d_r over K random unit directions; halves A/B (even/odd directions) give
# the noise correction: var(delta_true) ~ cov(delta_A, delta_B). Per layer and network:
#   X*    = cov(err, delta)^2 / (var(err) cov(delta_A, delta_B))   (noise-free share of the per-neuron error explained)
#   slope = cov(err, delta) / cov(delta_A, delta_B)                 (noise-free regression coefficient a in err ~ a delta)
#   rho   = std over directions of n d_r / |mean|                    (single-direction noise)
# Held-out merge: the noise-free slope of the training networks applied to the others, e - a delta_hat, with the noisy
# delta_hat (a lower bound on what an exact delta gives) and the noise-free prediction 1 - X*.
import glob, sys, numpy as np
R, T, H = sys.argv[1], sys.argv[2], sys.argv[3]; nets = [int(x) for x in sys.argv[4:]]
n = 1024; res = {}
for net in nets:
    fs = sorted(glob.glob(f"{R}/heatdef_{net}_*_h{H}.npz"), key=lambda f: int(f.rsplit("_", 2)[1]))
    if not fs:
        continue
    zs = [np.load(f) for f in fs]
    d = np.concatenate([z["d"] for z in zs]); idx = np.concatenate([z["idx"] for z in zs]); E0 = zs[0]["E0"]
    h = float(zs[0]["h"]); K = len(d)
    mt = np.load(f"{T}/truth_off{net}.npz")["m"].astype(np.float64)
    err = E0 - mt
    dA = 2 * n * d[idx % 2 == 0].mean(0); dB = 2 * n * d[idx % 2 == 1].mean(0); dl = 0.5 * (dA + dB)
    out = []
    for l in range(16):
        e = err[l] - err[l].mean(); a = dA[l] - dA[l].mean(); b = dB[l] - dB[l].mean(); m = dl[l] - dl[l].mean()
        cab = np.mean(a * b); ce = np.mean(e * m); ve = np.mean(e * e)
        xs = ce * ce / (ve * cab) if cab > 0 else np.nan
        rho = np.std(2 * n * d[:, l, :], axis=0).mean() / max(np.abs(dl[l]).mean(), 1e-30)
        out.append((np.mean(err[l] ** 2), xs, ce / cab if cab > 0 else np.nan, np.corrcoef(err[l], dl[l])[0, 1], rho,
                    np.sqrt(np.mean(dl[l] ** 2)), cab / np.mean(m * m)))
    res[net] = dict(out=np.array(out), err=err, dl=dl, K=K, h=h)
    print(f"net {net}: K = {K} directions, h = {h}; chain raw (this variant) {np.mean(err[-1] ** 2):.4e}")
print("\nlayer | MSE | X* (noise-free explained share) | slope a | raw corr | single-dir noise rho | rms delta | "
      "signal share of var(delta_hat)   [mean over networks]")
O = np.mean([res[k]["out"] for k in res], 0)
for l in range(16):
    print(f"{l:2d} | {O[l,0]:.2e} | {O[l,1]:+.3f} | {O[l,2]:+.4f} | {O[l,3]:+.3f} | {O[l,4]:.1f} | {O[l,5]:.2e} | {O[l,6]:.3f}")
ks = sorted(res)
if len(ks) >= 2:
    tr, te = ks[: len(ks) // 2], ks[len(ks) // 2:]
    for l in (12, 13, 14, 15):
        a = np.mean([res[k]["out"][l, 2] for k in tr])
        for k in te:
            e = res[k]["err"][l]; m = res[k]["dl"][l]
            print(f"held-out net {k} layer {l}: slope a = {a:+.4f} from {tr}; MSE {np.mean(e ** 2):.3e} -> "
                  f"{np.mean((e - a * (m - m.mean() * 0)) ** 2):.3e} with the noisy delta_hat "
                  f"({100 * (np.mean((e - a * m) ** 2) / np.mean(e ** 2) - 1):+.1f}%); noise-free X* here "
                  f"{res[k]['out'][l, 1]:.3f}")
