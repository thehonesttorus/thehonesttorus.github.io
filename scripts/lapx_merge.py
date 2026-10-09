"""Merge coordinate-Laplacian blocks (scripts/lap_exact.py) and evaluate the heat-defect correction.
  python scripts/lapx_merge.py RESULTS_DIR [kind=gauss] [nets=0,1] [h_ref=]"""
import sys, glob, numpy as np
d = sys.argv[1]; kw = dict(kind="gauss", nets="0,1")
for a in sys.argv[2:]:
    k, v = a.split("="); kw[k] = v
R = {}
for net in [int(x) for x in kw["nets"].split(",")]:
    fs = sorted(glob.glob(f"{d}/lapx_{kw['kind']}_{net}_c*.npz"))
    if not fs: continue
    zs = [np.load(f) for f in fs]; e0 = zs[0]["e0"]; truth = zs[0]["truth"]
    d2 = np.concatenate([z["d2"] for z in zs]); coords = np.concatenate([z["coords"] for z in zs]); n = len(e0)
    o = np.argsort(coords); d2 = d2[o]; coords = coords[o]; err = e0 - truth
    print(f"net {net} {kw['kind']}: {len(coords)}/{n} coordinates, raw MSE {np.mean(err**2):.4e}")
    if len(coords) < n: print("   incomplete; scaling the partial sum by n/len")
    lap = d2.sum(0) * (n / len(coords)); delta = e0 - lap; a = np.dot(delta, err) / np.dot(delta, delta)
    r2 = 1 - np.sum((err - a * delta) ** 2) / np.sum(err ** 2)
    print(f"   |lap| {np.sqrt(np.mean(lap**2)):.3e} |e0| {np.sqrt(np.mean(e0**2)):.3e} |delta| {np.sqrt(np.mean(delta**2)):.3e} |err| {np.sqrt(np.mean(err**2)):.3e}")
    print(f"   corr(err,delta) {np.corrcoef(err, delta)[0,1]:+.4f}  a {a:+.4f}  explained {r2:.3f}  merge MSE {np.mean((err - a*delta)**2):.4e}  midpoint MSE {np.mean(((e0+lap)/2-truth)**2):.3e}")
    # spectrum of the coordinate second derivatives (how traceless the Hessian is)
    hd = d2.std(0); print(f"   per-coordinate d2: rms over (coord, neuron) {np.sqrt(np.mean(d2**2)):.3e}; trace/n rms {np.sqrt(np.mean((lap/n)**2)):.3e}")
    R[net] = dict(err=err, delta=delta, a=a)
if len(R) >= 2:
    for i in R:
        for j in R:
            if i != j: print(f"   a[{i}]={R[i]['a']:+.4f} -> net {j}: merge MSE {np.mean((R[j]['err'] - R[i]['a']*R[j]['delta'])**2):.4e}")
