"""Merge Laplacian probes (scripts/lap_probe.py outputs) and evaluate the heat-defect correction.
  python scripts/lap_merge.py RESULTS_DIR [kind=gauss] [nets=0,1]
Reports, per net: raw MSE, corr(err, delta), fitted a, explained fraction, merge MSE vs number of probes; and the
cross-network calibration (a fitted on one net, applied to the other)."""
import sys, os, glob, numpy as np
d = sys.argv[1]; kw = dict(kind="gauss", nets="0,1")
for a in sys.argv[2:]:
    k, v = a.split("="); kw[k] = v
nets = [int(x) for x in kw["nets"].split(",")]; R = {}
for net in nets:
    fs = sorted(glob.glob(f"{d}/lap_{kw['kind']}_{net}_s*.npz"))
    if not fs: continue
    zs = [np.load(f) for f in fs]; e0 = zs[0]["e0"]; truth = zs[0]["truth"]
    d2 = np.concatenate([z["d2"] for z in zs]); seeds = np.concatenate([z["seeds"] for z in zs])
    o = np.argsort(seeds); d2 = d2[o]; R[net] = dict(e0=e0, truth=truth, d2=d2, err=e0 - truth)
    err = e0 - truth; print(f"net {net} {kw['kind']}: {len(d2)} probes, raw MSE {np.mean(err**2):.4e}")
    for S in [s for s in (2, 4, 8, 16, 32, 64, 128) if s <= len(d2)]:
        lap = d2[:S].mean(0); delta = e0 - lap; a = np.dot(delta, err) / np.dot(delta, delta)
        r2 = 1 - np.sum((err - a * delta) ** 2) / np.sum(err ** 2)
        # probe noise: spread of the S-probe Laplacian estimate (std of the mean over probes), relative to |delta|
        se = np.sqrt(np.mean(d2[:S].var(0, ddof=1) / S)) if S > 1 else np.nan
        print(f"   S={S:3d}: corr {np.corrcoef(err, delta)[0,1]:+.3f}  a {a:+.4f}  explained {r2:.3f}  merge MSE {np.mean((err - a*delta)**2):.4e}  "
              f"|delta| {np.sqrt(np.mean(delta**2)):.2e}  probe se {se:.2e}  midpoint MSE {np.mean(((e0+lap)/2 - truth)**2):.3e}")
if len(R) >= 2:
    print("cross-network calibration (a fitted on net i, applied to net j), all probes:")
    A = {}
    for net, r in R.items():
        lap = r["d2"].mean(0); r["delta"] = r["e0"] - lap; A[net] = np.dot(r["delta"], r["err"]) / np.dot(r["delta"], r["delta"])
    for i in R:
        for j in R:
            if i != j:
                rj = R[j]; print(f"   a[{i}]={A[i]:+.4f} -> net {j}: merge MSE {np.mean((rj['err'] - A[i]*rj['delta'])**2):.4e} (own a {A[j]:+.4f})")
